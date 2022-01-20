//#include "stdafx.h"
#include <chrono>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include <algorithm> // max_element, std::min
#include <utility> // pair
#include <math.h> // ceil

//using namespace std;

// for debugging only
#define OBSERVED_NODE 149

#define DEBUG 1

typedef std::vector<int> int_vector;
typedef std::vector<std::pair<int, int>> pair_vector;
typedef std::vector<std::vector<int>> interval_tree;
typedef std::vector<std::set<int>> set_vector;

// common graph data
struct graph_data_t {
    int biggest_chunk;
    int longest_neighbor_seq;
    int cluster_count;
} graph_summary;

__global__ void GeneratePairs(int* indexes, int* neighbors, int* vertex_start, int* pairs, int* vertex_length, int* cluster_sizes, int* cluster_starts, int vertex_count, int maxNum, int vertexOffset, int neighborOffset) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i >= maxNum) {
        return;
    }

    if(i == 0) {
        for (int j = 0; j < vertex_count; j++) {
            printf("%d: %d\n", j, indexes[j]);
        }
    }

    // search in an 'interval tree' to map onto correct vertex's combinations
    int index = vertex_count / 2;
    int array_start = 0;
    int array_size = vertex_count - 1;
    while (true) {
        // int normalized_index = vertex_length[vertexOffset] + indexes[index];
        // if (i == 13) {
        //     printf("Normalized index: %d\n", normalized_index);
        // }
        if (indexes[index] == i || (indexes[index] < i && indexes[index + 1] > i) || (indexes[index] < i && index == vertex_count - 1)) {
            for (int j = index; j < vertex_count - 2; j++){
                if (indexes[j] == indexes[j + 1]) {
                    index = j + 1;
                } else {
                    break;
                }
            }
            break;
        }

        if (indexes[index] < i) {
            array_start = index;
        } else {
            array_size /= 2;
        }
        // middle of half-array
        index = (2 * array_start + array_size) / 2; 
    }

    int pair_combination = i - indexes[index];

    // get neighbor indexes
    int edge_one = pairs[2 * pair_combination];
    int edge_two = pairs[2 * pair_combination + 1];

    // get neighbors
    int neighbor_one = neighbors[vertex_start[index] - neighborOffset + edge_one];
    int neighbor_two = neighbors[vertex_start[index] - neighborOffset + edge_two];

    int cluster_index = 3 * i;
    int cluster_size = vertex_length[neighbor_one] + vertex_length[neighbor_two] + vertex_length[vertexOffset + index] - 4;
    cluster_sizes[i] = cluster_size;

    //printf("Vertex: %d\n", index + vertexOffset);

    // if (index + vertexOffset == OBSERVED_NODE) {
    //     printf("%d: Length: %d, Neighbor one: %d, neighbor");
    // }

    // middle vertex and two neighbors
    cluster_starts[cluster_index] = vertexOffset + index;
    cluster_starts[cluster_index + 1] = neighbor_one;
    cluster_starts[cluster_index + 2] = neighbor_two;

    // debug print
    if ((neighbor_one == OBSERVED_NODE || neighbor_two == OBSERVED_NODE || vertexOffset + index == OBSERVED_NODE) && DEBUG) {
        printf("Cluster of size %d for nodes %d-%d-%d\n", cluster_size, neighbor_one, neighbor_two, vertexOffset + index);
        printf("Vertex start: %d, prev vertex start: %d, Neighbor offset: %d, edge_one: %d, edge_two: %d, pair_combination = %d\n", vertex_start[index], vertex_start[index-1], neighborOffset, edge_one, edge_two, pair_combination);
    }
}

__global__ void CountClusterExpectedForce(int* cluster_size, int* cluster_start, int* total_vertex_size, float* output, int maxNum, int vertex_count, int offset) {

    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= maxNum) {
        return;
    }

    int vertex = cluster_start[i];
    
    int size = cluster_size[i];
    int vertex_total = total_vertex_size[vertex];
    float normalized = (float) size/vertex_total;
    output[i] = -(__logf(normalized) * (normalized));
}

void read_graph(int_vector &vertex_start, int_vector &vertex_length, int_vector &neighbors, set_vector &neighbor_sets, std::string infilename, char delimiter = ' ', int ignore_weights = 0) 
{
	vertex_start.clear(); 
    vertex_length.clear();
    neighbors.clear();

	std::ifstream infile;
	infile.open(infilename);
	std::string temp;
	
	int last_node = -1;
	int node_count = 0;
    int vertex_count = 0;
    
    // start of vertex 0
    vertex_start.push_back(0);
    std::set<int> neighbor_set;

	while (getline(infile, temp, delimiter)) {
        int node = std::stoi(temp);
        node_count++;

		if (node != last_node) {
            if (node - last_node > 1) {
                std::cout << "ERROR " << node << std::endl;
            }
            if (vertex_count != 0) {
                vertex_length.push_back(vertex_count);
                vertex_start.push_back(node_count-1);
                neighbor_sets.push_back(neighbor_set);
            }
        	last_node = node;
            vertex_count = 0;
            neighbor_set = std::set<int>();
		}
		
        vertex_count++;

		getline(infile, temp);
        int neighbor = std::stoi(temp);
        if (neighbor == last_node) {
            // edge to itself, e.g. 1->1
            vertex_count--;
            node_count--;
            std::cout << "Skipping edge " << neighbor << " - " << neighbor << std::endl;
        } else {
            // normal case
            neighbors.push_back(neighbor); 
            neighbor_set.insert(neighbor);
        }
	}

    // length of the last vertex negihbors 
    vertex_length.push_back(vertex_count);
    neighbor_sets.push_back(neighbor_set);
}

void generate_pairs(int_vector &pairs, int_vector &pair_count, int highest_degree) {
    pairs.clear();
    pair_count.clear();

    // vertex degree 0
    pair_count.push_back(0);
    // vertex degree 1
    pair_count.push_back(0);

    int count = 0; 

    for (int i = 1; i < highest_degree; i++) {
        for (int j = 0; j < i; j++) {
            pairs.push_back(j);
            pairs.push_back(i);
            count++;
        }
        pair_count.push_back(count);
    }
}

void split_into_generating_chunks(int_vector &generating_chunks, int_vector &vertex_length, int_vector &pair_count, int_vector &chunk_size, int_vector &cluster_start, interval_tree &intervals, int blocks, int threads) {
    int cluster_count = 0;
    int size = 0;
    int chunk_start = 0;
    int current_vertices = 0;
    int limit = blocks * threads;
    int neighbors_in_chunk = 0;
    int biggest_chunk = 0;
    int longest_neighbor_seq = 0;
    int_vector indexes;

    for (int i = 0; i < vertex_length.size(); i++) {
        cluster_start.push_back(cluster_count);
        
        int vertex_combinations = pair_count[vertex_length[i]];
        cluster_count += vertex_combinations;
        if (size + vertex_combinations > limit) {
            if (current_vertices == 0) {
                std::invalid_argument("Insufficient size of blocks and threads for generating chunks. Required combined size at least " + std::to_string(vertex_combinations) + ".");
            }
            intervals.push_back(indexes);
            chunk_size.push_back(size);
            biggest_chunk = std::max(biggest_chunk, (int) indexes.size());
            longest_neighbor_seq = std::max(longest_neighbor_seq, neighbors_in_chunk);
            generating_chunks.push_back(chunk_start);
            chunk_start += current_vertices;
            generating_chunks.push_back(chunk_start);
            chunk_start++;
            current_vertices = 0;
            size = vertex_combinations;
            indexes = std::vector<int>(0);
            neighbors_in_chunk = vertex_length[i];
        } else {
            indexes.push_back(size);
            size += vertex_combinations;
            current_vertices++;
            neighbors_in_chunk += vertex_length[i];
        }
    }

    if (current_vertices != 0) {
        generating_chunks.push_back(chunk_start);
        generating_chunks.push_back(chunk_start + current_vertices - 1);
        intervals.push_back(indexes);
        chunk_size.push_back(size);
        biggest_chunk = std::max(biggest_chunk, (int) indexes.size());
        longest_neighbor_seq = std::max(longest_neighbor_seq, neighbors_in_chunk);
    }

    for (int i = 0; i < generating_chunks.size(); i += 2) {
        std::cout << generating_chunks[i] << "-" << generating_chunks[i+1] << " ";
    }
    std::cout << std::endl;

    graph_summary.biggest_chunk = biggest_chunk;
    graph_summary.longest_neighbor_seq = longest_neighbor_seq;
    graph_summary.cluster_count = cluster_count;
}

void check_error() {
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
        printf("Error!\nLast Error: %s\nDetails: %s\n", cudaGetErrorName(error), cudaGetErrorString(error));
        exit(error);
    }
}

int main(int argc, char* argv[]) { //takes a filename (es: fb_full) as input; print its ExF in result.txt 

	//cout << "This program determines the Expected Force of every node for each graph.\n Stores the results in 'FILENAME_results.txt'" << endl;
	
    int blocks = 4;
    int threads = 128;
    int streamCount = 2;

    int_vector vertex_start, vertex_length, neighbors, pairs, pair_count, generating_chunks, chunk_size, path_vertex_one, cluster_start;
    interval_tree intervals;
    set_vector neighbor_sets;

    std::string filename = argv[1];
    int ignore_weights = std::stoi(argv[2]);

    std::cout << "Evaluating file " << filename << std::endl;

    //int64_t duration;
    //int repetitions = 5;

    //for (int i = 0; i < repetitions; i++) {

        //reads graph
        read_graph(vertex_start, vertex_length, neighbors, neighbor_sets, filename, ' ', ignore_weights); //converts SNAP graph to sorted edgelist.

        int highest_degree = *std::max_element(vertex_length.begin(), vertex_length.end());

        generate_pairs(pairs, pair_count, highest_degree);

        std::cout << "Generating chunks" << std::endl;

        split_into_generating_chunks(generating_chunks, vertex_length, pair_count, chunk_size, cluster_start, intervals, blocks, threads);

        int biggest_chunk = graph_summary.biggest_chunk;
        int most_neighbors = graph_summary.longest_neighbor_seq;

        cudaStream_t streams[streamCount];
        std::vector<int*> index_pointers;
        std::vector<int*> vertex_start_pointers;
        std::vector<int*> neighbor_pointers;
        std::vector<int*> cluster_size_pointers;
        std::vector<int*> host_cluster_size_pointers;
        std::vector<int*> cluster_start_pointers;
        std::vector<int*> host_cluster_start_pointers;
        for (int i = 0; i < streamCount; i++) {
            cudaStreamCreate(&streams[i]);
            std::cout << "stream created" << std::endl;
            check_error();

            int* index_ptr;
            cudaMalloc((void**)&index_ptr, sizeof(int) * biggest_chunk);
            index_pointers.push_back(index_ptr);
            
            int* vertex_start_ptr;
            std::cout << "Allocating " << sizeof(int) * biggest_chunk << "B of data for vertex_start_ptr" << std::endl;
            cudaMalloc((void**)&vertex_start_ptr, sizeof(int) * biggest_chunk);
            vertex_start_pointers.push_back(vertex_start_ptr);
            std::cout << "Allocated for vertex start: " << sizeof(int) * biggest_chunk << std::endl;

            int* neighbor_ptr;
            cudaMalloc((void**)&neighbor_ptr, sizeof(int) * most_neighbors);
            neighbor_pointers.push_back(neighbor_ptr);
            std::cout << "Most neighbors: " << most_neighbors << std::endl;

            int* cluster_size_ptr;
            cudaMalloc((void**)&cluster_size_ptr, sizeof(int) * graph_summary.cluster_count);
            cluster_size_pointers.push_back(cluster_size_ptr);

            int* host_cluster_size_ptr = (int*) malloc(sizeof(int) * graph_summary.cluster_count);
            host_cluster_size_pointers.push_back(host_cluster_size_ptr);

            int* cluster_start_ptr;
            cudaMalloc((void**)&cluster_start_ptr, sizeof(int) * graph_summary.cluster_count * 3);
            cluster_start_pointers.push_back(cluster_start_ptr);

            int* host_cluster_start_ptr = (int*) malloc(sizeof(int) * graph_summary.cluster_count * 3);
            host_cluster_start_pointers.push_back(host_cluster_start_ptr);
        }
        std::cout << "Pointer vectors allocated" << std::endl;
        check_error();

        int* pairs_ptr;
        cudaMalloc((void**)&pairs_ptr, sizeof(int) * pairs.size());
        cudaMemcpy(pairs_ptr, pairs.data(), sizeof(int) * pairs.size(), cudaMemcpyHostToDevice);
        std::cout << "Pairs pointer allocated" << std::endl;
        check_error();

        int* length_ptr;
        cudaMalloc((void**)&length_ptr, sizeof(int) * vertex_length.size());
        cudaMemcpy(length_ptr, vertex_length.data(), sizeof(int) * vertex_length.size(), cudaMemcpyHostToDevice);
        std::cout << "Length pointer allocated" << std::endl;
        check_error();

        std::vector<int> cluster_sizes;
        std::vector<int> start_vertex;
        std::vector<int> total_cluster_size(vertex_length.size(), 0);
        
        std::cout << "Number of chunks: " << generating_chunks.size() << std::endl;

        for (int index = 0; index < generating_chunks.size(); index += 2 * streamCount) {
            std::cout << "Computing chunk " << index / 2 + 1 << " of " << generating_chunks.size() / 2 << std::endl;
            int streamsUsed = std::min((int) (generating_chunks.size() / 2) - index/2, streamCount);
            std::cout << "Streams used: " << streamsUsed << std::endl;
            for (int i = 0; i < streamsUsed; i++) {
                int chunk_index = index + 2 * i;
                int chunk_start = generating_chunks[chunk_index];
                int chunk_end = generating_chunks[chunk_index + 1];
                int number_of_clusters = cluster_start[chunk_end] + pair_count[vertex_length[chunk_end]] - cluster_start[chunk_start];
                std::cout << "chunk start " << chunk_start << std::endl;
                std::cout << "chunk end " << chunk_end << std::endl;
                std::cout << "last chunk start " << vertex_start[chunk_end] << std::endl;
                std::cout << "last neighbor " << vertex_start[chunk_end] + vertex_length[chunk_end] - 1 << std::endl;
                int neighbor_size = vertex_start[chunk_end] + vertex_length[chunk_end] - vertex_start[chunk_start];
                std:: cout << "Neighbor size: " << neighbor_size << std::endl;

                cudaMemcpyAsync(index_pointers[i], intervals[index/2].data(), sizeof(int) * intervals[index/2].size(), cudaMemcpyHostToDevice);
                std::cout << "Indexes copied" << std::endl;
                check_error();

                std::cout << "Copying " << sizeof(int) * (chunk_end - chunk_start + 1) << "B of data" << std::endl;
                std::cout << "Start of last vertex: " << vertex_start[chunk_end] << std::endl;
                cudaMemcpy(vertex_start_pointers[i], vertex_start.data() + chunk_start, sizeof(int) * (chunk_end - chunk_start + 1), cudaMemcpyHostToDevice);
                std::cout << "Vertex starts copied" << std::endl;
                check_error();

                std::cout << "Copying vertex starts from " << vertex_start[chunk_start] << " with length of " << neighbor_size << std::endl;
                cudaMemcpy(neighbor_pointers[i], neighbors.data() + vertex_start[chunk_start], sizeof(int) * neighbor_size, cudaMemcpyHostToDevice);
                std::cout << "Neighbors copied" << std::endl;
                check_error();

                GeneratePairs<<<blocks, threads, 0, streams[i]>>>(index_pointers[i], neighbor_pointers[i], vertex_start_pointers[i], pairs_ptr, length_ptr, cluster_size_pointers[i], cluster_start_pointers[i], chunk_end - chunk_start + 1, number_of_clusters, chunk_start, vertex_start[chunk_start]);
                cudaDeviceSynchronize();
                std::cout << "Generation finished" << std::endl;
                check_error();

                std::cout << "Copieeed clusters: " << number_of_clusters << std::endl;
                cudaMemcpyAsync(host_cluster_size_pointers[i], cluster_size_pointers[i], sizeof(int) * number_of_clusters, cudaMemcpyDeviceToHost);
                std::cout << "host cluster sizes copied" << std::endl;
                check_error();

                cudaMemcpyAsync(host_cluster_start_pointers[i], cluster_start_pointers[i], sizeof(int) * number_of_clusters * 3, cudaMemcpyDeviceToHost);
                std::cout << "host cluster starts copied" << std::endl;
                check_error();
            }

            std::cout << "Pairs Generated" << std::endl;
            check_error();

            cudaDeviceSynchronize();

            for (int i = 0; i < streamsUsed; i++) {
                int* path_ptr = host_cluster_start_pointers[i];
                int chunk_index = 2 * (index + i);
                int chunk_start = generating_chunks[chunk_index];
                int chunk_end = generating_chunks[chunk_index + 1];
                int number_of_clusters = cluster_start[chunk_end] + pair_count[vertex_length[chunk_end]] - cluster_start[chunk_start];

                for (int j = 0; j < number_of_clusters; j++) {
                    int cluster_size = host_cluster_size_pointers[i][j];
                    int source_vertex = path_ptr[3 * j];
                    int neighbor_one = path_ptr[3 * j + 1];
                    int neighbor_two = path_ptr[3 * j + 2];

                    // if two neighbors form a triangle with cluster-size edges
                    // contains() cannot be used -> C++20-only
                    if (neighbor_sets[neighbor_one].find(neighbor_two) != neighbor_sets[neighbor_one].end()) {
                        cluster_size -= 2;
                    }

                    if (source_vertex == OBSERVED_NODE && DEBUG) {
                        std::cout << "Cluster " << OBSERVED_NODE << ", size: " << cluster_size << std::endl;
                    }

                    //std::cout << "Processing cluster starting at " << source_vertex << " with size " << cluster_size << std::endl;
                    cluster_sizes.push_back(cluster_size);
                    cluster_sizes.push_back(cluster_size);
                    cluster_sizes.push_back(cluster_size);
                    cluster_sizes.push_back(cluster_size);

                    // push twice - two combinations of S->A S->B and S->B S->A (two clusters but with same edges) 
                    start_vertex.push_back(source_vertex);
                    start_vertex.push_back(source_vertex);
                    start_vertex.push_back(neighbor_one);
                    start_vertex.push_back(neighbor_two);
                    total_cluster_size[source_vertex] += cluster_size;
                    total_cluster_size[source_vertex] += cluster_size;
                    total_cluster_size[neighbor_one] += cluster_size;
                    total_cluster_size[neighbor_two] += cluster_size;
                }
            }
        }

        std::cout << "Total cluster sizes: " << std::endl;
        for (int i = 0; i < total_cluster_size.size(); i++) {
            std::cout << i << ": " << total_cluster_size[i] << std::endl;
        }

        for (int i = 0; i < start_vertex.size(); i++) {
            if (start_vertex[i] == OBSERVED_NODE && DEBUG) {
                std::cout << "i: " << i << ", Cluster size:" << cluster_sizes[i] << " total: " << total_cluster_size[OBSERVED_NODE] << std::endl;
            }
        }

        std::cout << "end" << std::endl;

        for (int i = 0; i < streamCount; i++) {
            cudaFree(index_pointers[i]);
            cudaFree(vertex_start_pointers[i]);
            cudaFree(neighbor_pointers[i]);
            cudaFree(cluster_size_pointers[i]);
            cudaFree(cluster_start_pointers[i]);
            free(host_cluster_size_pointers[i]);
            free(host_cluster_start_pointers[i]);
        }     
        cudaFree(length_ptr);
        cudaFree(pairs_ptr);

        std::vector<int*> input_sizes;
        std::vector<int*> input_vertices;
        std::vector<float*> outputs;
        int array_size = blocks * threads;
        for (int i = 0; i < streamCount; i++) {
            int *size_ptr, *vertex_ptr; 
            float *devOut;
            cudaMalloc((void**)&size_ptr, sizeof(int) * array_size);
            cudaMalloc((void**)&vertex_ptr, sizeof(int) * array_size);
            cudaMalloc((void**)&devOut, sizeof(float) * array_size);
            input_sizes.push_back(size_ptr);
            input_vertices.push_back(vertex_ptr);
            outputs.push_back(devOut);
        }

        int *total_size_ptr;
        cudaMalloc((void**)&total_size_ptr, sizeof(int) * total_cluster_size.size());

        cudaMemcpyAsync(total_size_ptr, total_cluster_size.data(), sizeof(int) * total_cluster_size.size(), cudaMemcpyHostToDevice);

        // std::ceil is stupid
        int ceil_share = (int) ((graph_summary.cluster_count * 4 / array_size) + ((graph_summary.cluster_count * 4 % array_size) != 0));
        std::cout << "Clusters total: " << graph_summary.cluster_count * 4 << ", array_size: " << array_size << std::endl;
        std::cout << "Ratio: " << ((graph_summary.cluster_count * 4) / array_size) << ", ceil: " << ceil_share << std::endl;
        int chunks = std::max(1, ceil_share);

        std::cout << "Computing normalized force in " << chunks << " chunks" << std::endl;

        std::vector<float> normalized_sizes(cluster_sizes.size(), 0);

        for (int index = 0; index < chunks; index += streamCount) {
            int streamsUsed = std::min(streamCount, chunks - index);
            for (int i = 0; i < streamsUsed; i++) {
                int element_count = std::min(array_size, 4 * graph_summary.cluster_count - i * array_size);
                cudaMemcpyAsync(input_sizes[i], cluster_sizes.data() + i * array_size, sizeof(int) * element_count, cudaMemcpyHostToDevice);
                cudaMemcpyAsync(input_vertices[i], start_vertex.data() + i * array_size, sizeof(int) * element_count, cudaMemcpyHostToDevice);
                CountClusterExpectedForce<<<blocks, threads, 0, streams[i]>>>(input_sizes[i], input_vertices[i], total_size_ptr, outputs[i], element_count, vertex_length.size(), array_size * i);
                cudaMemcpyAsync(normalized_sizes.data() + i * array_size, outputs[i], sizeof(float) * element_count, cudaMemcpyDeviceToHost);
            }

            cudaDeviceSynchronize();

            check_error();
        }

        std::cout << "Normalized sizes counted, freeing cuda arrays: " << input_sizes.size() << std::endl;

        for (int i = 0; i < streamCount; i++) {
            cudaFree(input_sizes[i]);
            check_error();

            cudaFree(input_vertices[i]);
            check_error();

            cudaFree(outputs[i]);
            check_error();
        }    

        cudaFree(total_size_ptr);

        // TODO: deallocate streams

        std::cout << "CUDA pointers deallocated" << std::endl;

        std::vector<float> results(vertex_length.size(), 0);

        std::cout << "Starting summation of normalized edges\n" << std::endl;

        for (int i = 0; i < normalized_sizes.size() - 1; i++) {
            results[start_vertex[i]] += normalized_sizes[i];
        }

        for (int i = 0; i < results.size(); i++) {
            std::cout << i << ": " << results[i] << std::endl;
        }

	return 0;
}
