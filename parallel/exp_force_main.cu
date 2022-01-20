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

// for debugging only
#define OBSERVED_NODE -1
#define DEBUG 0

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

/**
 * Kernel for generating all paths of length 2 in the graph and calculating the cluster sizes.
 */
__global__ void GeneratePairs(int* indexes, int* neighbors, int* vertex_start, int* pairs, int* vertex_length, int* cluster_sizes, int* cluster_starts, int vertex_count, int maxNum, int vertexOffset, int neighborOffset) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i >= maxNum) {
        return;
    }

    // search in an 'interval tree' to map onto correct vertex's combinations
    int index = 0;
    int array_start = 0;
    int array_size = vertex_count - 1;
    while (array_start <= array_size) {
        index = (array_start + array_size) / 2;
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
            array_start = index + 1;
        } else {
            array_size = index - 1;
        }
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

    // middle vertex and two neighbors
    cluster_starts[cluster_index] = vertexOffset + index;
    cluster_starts[cluster_index + 1] = neighbor_one;
    cluster_starts[cluster_index + 2] = neighbor_two;
}

/**
 * Calculates the expected force for the cluster.
 */
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

/**
 * Reads the graph.
 */
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

/**
 * Generates pair combinations in 2 values i-j, where i is i-th neighbor of source vertex and j is the source vertex's j-th neighbor.
 */
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

/**
 * Splits all vertices in graph into chunks that can be computed at once. Also calculates several other properties of the graph.
 */
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
        if (size + vertex_combinations > limit || vertex_combinations > limit) {
            if (current_vertices == 0) {
                std::invalid_argument("Insufficient size of blocks and threads for generating chunks. Required combined size at least " + std::to_string(vertex_combinations) + ".");
            }
            intervals.push_back(indexes);
            chunk_size.push_back(size);
            biggest_chunk = std::max(biggest_chunk, (int) indexes.size());
            longest_neighbor_seq = std::max(longest_neighbor_seq, neighbors_in_chunk);
            generating_chunks.push_back(chunk_start);
            chunk_start += current_vertices;
            generating_chunks.push_back(chunk_start - 1);
            current_vertices = 1;
            size = vertex_combinations;
            indexes = std::vector<int>(1, 0);
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

    graph_summary.biggest_chunk = biggest_chunk;
    graph_summary.longest_neighbor_seq = longest_neighbor_seq;
    graph_summary.cluster_count = cluster_count;
}

/**
 * Checks for cuda errors, logs if there is any error and finishes execution.
 */
void check_error() {
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
        printf("Error!\nLast Error: %s\nDetails: %s\n", cudaGetErrorName(error), cudaGetErrorString(error));
        exit(error);
    }
}

/**
 * Entry function, parameters: 
 * - filename: name fo the file whence the graph shall be read
 * - blocks - nubmer of blocks to be used
 * - threads - number of threads per block to be used, max. 1024 (GPU limitation)
 * - streams - 1..n
 */
int main(int argc, char* argv[]) { //takes a filename (es: fb_full) as input; print its ExF in result.txt 

	//cout << "This program determines the Expected Force of every node for each graph.\n Stores the results in 'FILENAME_results.txt'" << endl;
	
    if(argc < 5) {
        std::cout << "Insufficient number of arguments: " << argc << std::endl;
        exit(3);
    }

    std::string filename = argv[1];
    int blocks = atoi(argv[2]);
    int threads = atoi(argv[3]);
    int streamCount = atoi(argv[4]);

    int_vector vertex_start, vertex_length, neighbors, pairs, pair_count, generating_chunks, chunk_size, path_vertex_one, cluster_start;
    interval_tree intervals;
    set_vector neighbor_sets;

    //int ignore_weights = std::stoi(argv[2]);

    std::cout << "Evaluating file " << filename << std::endl;

    int64_t duration;
    int repetitions = 5;

    //reads graph
    read_graph(vertex_start, vertex_length, neighbors, neighbor_sets, filename, ' ', 1); //converts graph to a v-graph-like structure

    for (int i = 0; i < repetitions; i++) {

        auto start = std::chrono::high_resolution_clock::now();

        int highest_degree = *std::max_element(vertex_length.begin(), vertex_length.end());

        generate_pairs(pairs, pair_count, highest_degree);

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
            check_error();

            int* index_ptr;
            cudaMalloc((void**)&index_ptr, sizeof(int) * biggest_chunk);
            index_pointers.push_back(index_ptr);
            
            int* vertex_start_ptr;
            cudaMalloc((void**)&vertex_start_ptr, sizeof(int) * biggest_chunk);
            vertex_start_pointers.push_back(vertex_start_ptr);

            int* neighbor_ptr;
            cudaMalloc((void**)&neighbor_ptr, sizeof(int) * most_neighbors);
            neighbor_pointers.push_back(neighbor_ptr);

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
        check_error();

        int* pairs_ptr;
        cudaMalloc((void**)&pairs_ptr, sizeof(int) * pairs.size());
        cudaMemcpy(pairs_ptr, pairs.data(), sizeof(int) * pairs.size(), cudaMemcpyHostToDevice);
        check_error();

        int* length_ptr;
        cudaMalloc((void**)&length_ptr, sizeof(int) * vertex_length.size());
        cudaMemcpy(length_ptr, vertex_length.data(), sizeof(int) * vertex_length.size(), cudaMemcpyHostToDevice);
        check_error();

        std::vector<int> cluster_sizes;
        std::vector<int> start_vertex;
        std::vector<int> total_cluster_size(vertex_length.size(), 0);

        for (int index = 0; index < generating_chunks.size(); index += 2 * streamCount) {
            int streamsUsed = std::min((int) (generating_chunks.size() / 2) - index/2, streamCount);
            for (int i = 0; i < streamsUsed; i++) {
                int chunk_index = index + 2 * i;
                int chunk_start = generating_chunks[chunk_index];
                int chunk_end = generating_chunks[chunk_index + 1];
                int number_of_clusters = cluster_start[chunk_end] + pair_count[vertex_length[chunk_end]] - cluster_start[chunk_start];
                int neighbor_size = vertex_start[chunk_end] + vertex_length[chunk_end] - vertex_start[chunk_start];

                cudaMemcpyAsync(index_pointers[i], intervals[index/2 + i].data(), sizeof(int) * intervals[index/2 + i].size(), cudaMemcpyHostToDevice);
                check_error();

                cudaMemcpy(vertex_start_pointers[i], vertex_start.data() + chunk_start, sizeof(int) * (chunk_end - chunk_start + 1), cudaMemcpyHostToDevice);
                check_error();

                cudaMemcpy(neighbor_pointers[i], neighbors.data() + vertex_start[chunk_start], sizeof(int) * neighbor_size, cudaMemcpyHostToDevice);
                check_error();

                GeneratePairs<<<blocks, threads, 0, streams[i]>>>(index_pointers[i], neighbor_pointers[i], vertex_start_pointers[i], pairs_ptr, length_ptr, cluster_size_pointers[i], cluster_start_pointers[i], chunk_end - chunk_start + 1, number_of_clusters, chunk_start, vertex_start[chunk_start]);
                cudaDeviceSynchronize();
                check_error();

                cudaMemcpy(host_cluster_size_pointers[i], cluster_size_pointers[i], sizeof(int) * number_of_clusters, cudaMemcpyDeviceToHost);
                check_error();

                cudaMemcpy(host_cluster_start_pointers[i], cluster_start_pointers[i], sizeof(int) * number_of_clusters * 3, cudaMemcpyDeviceToHost);
                check_error();
            }
            check_error();

            cudaDeviceSynchronize();

            for (int i = 0; i < streamsUsed; i++) {
                int* path_ptr = host_cluster_start_pointers[i];
                int chunk_index = index + 2 * i;
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
        int chunks = std::max(1, ceil_share);

        std::vector<float> normalized_sizes(cluster_sizes.size(), 0);

        for (int index = 0; index < chunks; index += streamCount) {
            int streamsUsed = std::min(streamCount, chunks - index);
            for (int i = 0; i < streamsUsed; i++) {
                int element_count = std::min(array_size, 4 * graph_summary.cluster_count - (index + i) * array_size);
                cudaMemcpyAsync(input_sizes[i], cluster_sizes.data() + (index + i) * array_size, sizeof(int) * element_count, cudaMemcpyHostToDevice);
                cudaMemcpyAsync(input_vertices[i], start_vertex.data() + (index + i) * array_size, sizeof(int) * element_count, cudaMemcpyHostToDevice);
                CountClusterExpectedForce<<<blocks, threads, 0, streams[i]>>>(input_sizes[i], input_vertices[i], total_size_ptr, outputs[i], element_count, vertex_length.size(), array_size * i);
                cudaMemcpyAsync(normalized_sizes.data() + (index + i) * array_size, outputs[i], sizeof(float) * element_count, cudaMemcpyDeviceToHost);
            }

            cudaDeviceSynchronize();
            check_error();
        }

        for (int i = 0; i < streamCount; i++) {
            cudaFree(input_sizes[i]);
            check_error();

            cudaFree(input_vertices[i]);
            check_error();

            cudaFree(outputs[i]);
            check_error();

            cudaStreamDestroy(streams[i]);
            check_error();
        }    
        cudaFree(total_size_ptr);

        std::vector<float> results(vertex_length.size(), 0);

        for (int i = 0; i < normalized_sizes.size(); i++) {
            results[start_vertex[i]] += normalized_sizes[i];
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        duration += time;
    }

    std::cout << "Time in ms (" << blocks << ", " << threads << ", " << streamCount"):" << std::endl << duration/repetitions << std::endl;

    std::ofstream outfile;
    outfile.open("results_" + filename);

    for (int i = 0; i < results.size(); i++) {
        outfile << i << "  " << results[i] << std::endl;
    }

	return 0;
}
