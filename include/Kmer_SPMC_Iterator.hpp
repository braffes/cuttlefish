#ifndef KMER_SPMC_ITERATOR_HPP
#define KMER_SPMC_ITERATOR_HPP

#include "Kmer.hpp"
#include "Kmer_Container.hpp"
#include "kmc_api/kmc_file.h"

#include <atomic>
#include <cstdint>
#include <cstddef>
#include <memory>
#include <vector>
#include <string>
#include <thread>
#include <iostream>
#include <chrono>

// Data required by the consumers to correctly parse raw binary k-mers.
struct alignas(64)  // Use 64 for cache alignment if L1_CACHE_LINE_SIZE not defined
    Consumer_Data
{
    uint8_t* suff_buf{nullptr};
    uint64_t kmers_available;
    uint64_t kmers_parsed;
    std::vector<std::pair<uint64_t, uint64_t>> pref_buf;
    std::vector<std::pair<uint64_t, uint64_t>>::iterator pref_it;
};

template <uint16_t k>
class Kmer_SPMC_Iterator
{
private:
    const Kmer_Container<k>* const kmer_container;
    CKMC_DB kmer_database;
    const uint64_t kmer_count;
    const size_t consumer_count;
    uint64_t kmers_read;

    std::unique_ptr<std::thread> reader{nullptr};
    
    // Optimized buffer sizes
    static constexpr size_t BUF_SZ_PER_CONSUMER = (1 << 24);  // 64 MB for better EBS performance
    static constexpr size_t PREFETCH_BUFFER_COUNT = 20;        // Double buffering

    std::vector<Consumer_Data> consumer;
    std::vector<std::vector<uint8_t>> suffix_buffers;
    
    enum class Task_Status: uint8_t { pending, available, no_more };
    std::atomic<Task_Status>* task_status{nullptr};

    void read_raw_kmers();
    size_t get_idle_consumer() const;
    void open_kmer_database(const std::string& db_path);
    void close_kmer_database();

public:
    // Keep the original interface exactly the same
    Kmer_SPMC_Iterator(const Kmer_Container<k>* kmer_container, 
                      size_t consumer_count, 
                      bool at_begin = true, 
                      bool at_end = false);

    Kmer_SPMC_Iterator(const Kmer_SPMC_Iterator& other);
    ~Kmer_SPMC_Iterator();

    Kmer_SPMC_Iterator& operator=(const Kmer_SPMC_Iterator& rhs) = delete;

    bool value_at(size_t consumer_id, Kmer<k>& kmer);
    bool operator==(const Kmer_SPMC_Iterator& rhs) const;
    bool operator!=(const Kmer_SPMC_Iterator& rhs) const;
    void launch_production();
    void seize_production();
    bool launched() const;
    bool tasks_expected(size_t consumer_id) const;
    bool task_available(size_t consumer_id) const;
    std::size_t memory() const;
    static std::size_t memory(std::size_t consumer_count);

    // Dummy methods
    const Kmer_SPMC_Iterator& operator++() { return *this; }
    Kmer<k> operator*() { return Kmer<k>(); }
};

// Implementation with original class name but optimized internals

template <uint16_t k>
inline Kmer_SPMC_Iterator<k>::Kmer_SPMC_Iterator(
    const Kmer_Container<k>* const kmer_container, 
    const size_t consumer_count, 
    const bool at_begin, 
    const bool at_end)
    : kmer_container(kmer_container),
      kmer_count{kmer_container->size()},
      consumer_count{consumer_count},
      kmers_read{at_end ? kmer_count : 0}
{
    if(!(at_begin ^ at_end))
    {
        std::cerr << "Invalid position provided for SPMC k-mer iterator. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}

template <uint16_t k>
inline Kmer_SPMC_Iterator<k>::Kmer_SPMC_Iterator(const Kmer_SPMC_Iterator& other)
    : kmer_container(other.kmer_container),
      kmer_count{other.kmer_count},
      consumer_count{other.consumer_count},
      kmers_read{other.kmers_read}
{}

template <uint16_t k>
inline Kmer_SPMC_Iterator<k>::~Kmer_SPMC_Iterator()
{
    if(task_status != nullptr)
    {
        delete[] task_status;
        std::cerr << "\nCompleted a pass over the k-mer database.\n";
    }
    
    if(reader && reader->joinable())
    {
        reader->join();
    }
}

template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::open_kmer_database(const std::string& db_path)
{
    if(!kmer_database.open_for_cuttlefish_listing(db_path))
    {
        std::cerr << "Error opening k-mer database with prefix " << db_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}

template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::close_kmer_database()
{
    if(!kmer_database.Close())
    {
        std::cerr << "Error closing k-mer database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}

template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::launch_production()
{
    if(launched())
        return;

    // Initialize with optimized buffers
    task_status = new std::atomic<Task_Status>[consumer_count];
    consumer.resize(consumer_count);
    suffix_buffers.resize(consumer_count * PREFETCH_BUFFER_COUNT);

    for(size_t id = 0; id < consumer_count; ++id)
    {
        // Allocate multiple buffers for prefetching
        for(size_t buf_idx = 0; buf_idx < PREFETCH_BUFFER_COUNT; ++buf_idx)
        {
            suffix_buffers[id * PREFETCH_BUFFER_COUNT + buf_idx].resize(BUF_SZ_PER_CONSUMER);
        }
        
        consumer[id].suff_buf = suffix_buffers[id * PREFETCH_BUFFER_COUNT].data();
        consumer[id].kmers_available = 0;
        consumer[id].kmers_parsed = 0;
        consumer[id].pref_buf.clear();
        consumer[id].pref_it = consumer[id].pref_buf.begin();
        task_status[id] = Task_Status::pending;
    }

    open_kmer_database(kmer_container->container_location());

    reader.reset(new std::thread([this]() { read_raw_kmers(); }));
}

template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::launched() const
{
    return reader != nullptr;
}

template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::read_raw_kmers()
{
    std::vector<size_t> current_buffer_index(consumer_count, 0);
    
    while(!kmer_database.Eof())
    {
        bool found_idle_consumer = false;
        
        for(size_t id = 0; id < consumer_count; ++id)
        {
            if(task_status[id] == Task_Status::pending)
            {
                found_idle_consumer = true;
                Consumer_Data& consumer_state = consumer[id];
                
                // Switch to next buffer
                size_t next_buf_idx = (current_buffer_index[id] + 1) % PREFETCH_BUFFER_COUNT;
                consumer_state.suff_buf = suffix_buffers[id * PREFETCH_BUFFER_COUNT + next_buf_idx].data();
                
                consumer_state.kmers_available = kmer_database.read_raw_suffixes(
                    consumer_state.suff_buf, 
                    consumer_state.pref_buf, 
                    BUF_SZ_PER_CONSUMER
                );
                
                if(!consumer_state.kmers_available && !kmer_database.Eof())
                {
                    std::cerr << "Error reading the suffix file. Aborting.\n";
                    std::exit(EXIT_FAILURE);
                }

                kmers_read += consumer_state.kmers_available;
                consumer_state.kmers_parsed = 0;
                consumer_state.pref_it = consumer_state.pref_buf.begin();
                task_status[id] = Task_Status::available;
                current_buffer_index[id] = next_buf_idx;
                
                if(kmer_database.Eof()) break;
            }
        }
        
        // Prevent busy waiting
        if(!found_idle_consumer && !kmer_database.Eof())
        {
            std::this_thread::sleep_for(std::chrono::microseconds(100));
        }
        
        if(kmer_database.Eof()) break;
    }
    
    // Wait for consumers to finish current work
    for(size_t id = 0; id < consumer_count; ++id)
    {
        while(task_status[id] != Task_Status::pending)
        {
            std::this_thread::sleep_for(std::chrono::microseconds(10));
        }
    }
}

template <uint16_t k>
inline size_t Kmer_SPMC_Iterator<k>::get_idle_consumer() const
{
    size_t id = 0;
    while(task_status[id] != Task_Status::pending)
    {
        id = (id + 1) % consumer_count;
    }
    return id;
}

template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::seize_production()
{
    if(reader && reader->joinable())
    {
        reader->join();
    }

    for(size_t id = 0; id < consumer_count; ++id)
    {
        task_status[id] = Task_Status::no_more;
    }

    close_kmer_database();
}

template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::value_at(size_t consumer_id, Kmer<k>& kmer)
{
    if(!task_available(consumer_id))
        return false;

    auto& ts = consumer[consumer_id];
    if(ts.kmers_parsed == ts.kmers_available)
    {
        task_status[consumer_id] = Task_Status::pending;
        return false;
    }

    kmer_database.parse_kmer_buf<k>(ts.pref_it, ts.suff_buf, 
                                  ts.kmers_parsed * kmer_database.suff_record_size(), 
                                  kmer);
    ts.kmers_parsed++;

    return true;
}

template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::operator==(const Kmer_SPMC_Iterator& rhs) const
{
    return kmer_container == rhs.kmer_container && kmers_read == rhs.kmers_read;
}

template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::operator!=(const Kmer_SPMC_Iterator& rhs) const
{
    return !operator==(rhs);
}

template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::tasks_expected(size_t consumer_id) const
{
    return task_status[consumer_id] != Task_Status::no_more;
}

template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::task_available(size_t consumer_id) const
{
    return task_status[consumer_id] == Task_Status::available;
}

template <uint16_t k>
inline std::size_t Kmer_SPMC_Iterator<k>::memory() const
{
    std::size_t total = 0;
    total += consumer_count * sizeof(std::atomic<Task_Status>);
    total += consumer_count * sizeof(Consumer_Data);
    total += consumer_count * PREFETCH_BUFFER_COUNT * BUF_SZ_PER_CONSUMER;
    
    for(const auto& c : consumer)
    {
        total += c.pref_buf.capacity() * sizeof(std::pair<uint64_t, uint64_t>);
    }
    
    return total + CKMC_DB::pref_buf_memory();
}

template <uint16_t k>
inline std::size_t Kmer_SPMC_Iterator<k>::memory(std::size_t consumer_count)
{
    std::size_t total = 0;
    total += consumer_count * sizeof(std::atomic<Task_Status>);
    total += consumer_count * sizeof(Consumer_Data);
    total += consumer_count * PREFETCH_BUFFER_COUNT * BUF_SZ_PER_CONSUMER;
    total += consumer_count * (1 << 20) * sizeof(std::pair<uint64_t, uint64_t>);
    
    return total + CKMC_DB::pref_buf_memory();
}

#endif // KMER_SPMC_ITERATOR_HPP
