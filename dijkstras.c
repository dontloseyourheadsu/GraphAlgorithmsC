#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define MAX_VERTICES 100

struct MinHeapNode
{
    int v; // Vertex number
    int dist; // Distance of this vertex from source
};

struct MinHeap
{
    int size; // Number of heap nodes present currently
    int capacity; // Capacity of min heap
    struct MinHeapNode **array; // This is needed for decreaseKey()
    int *pos; // Position of nodes in the min heap
};

// A utility function to create a new Min Heap Node
struct MinHeapNode *newMinHeapNode(int v, int dist)
{
    struct MinHeapNode *minHeapNode = (struct MinHeapNode *)malloc(sizeof(struct MinHeapNode)); // Allocate memory for the node
    minHeapNode->v = v; // Set the vertex number
    minHeapNode->dist = dist; // Set the distance
    return minHeapNode; // Return the node
}

// A utility function to create a new adjacency matrix
int **createGraph(int vertices)
{
    int **graph = (int **)malloc(vertices * sizeof(int *));
    for (int i = 0; i < vertices; ++i)
    {
        graph[i] = (int *)malloc(vertices * sizeof(int));
        for (int j = 0; j < vertices; ++j)
        {
            graph[i][j] = (i == j) ? 0 : INT_MAX;
        }
    }
    return graph;
}

// A utility function to create a Min Heap
struct MinHeap *createMinHeap(int capacity)
{
    struct MinHeap *minHeap = (struct MinHeap *)malloc(sizeof(struct MinHeap)); // Allocate memory for the heap
    minHeap->pos = (int *)malloc(capacity * sizeof(int)); // Allocate memory for the position array
    minHeap->size = 0; // Initialize size to 0
    minHeap->capacity = capacity; // Set capacity
    minHeap->array = (struct MinHeapNode **)malloc(capacity * sizeof(struct MinHeapNode *)); // Allocate memory for the array
    return minHeap; // Return the heap
}

// A utility function to swap two nodes of min heap. Needed for min heapify
void swapMinHeapNode(struct MinHeapNode **a, struct MinHeapNode **b)
{
    struct MinHeapNode *t = *a; // Store the first node in a temporary variable
    *a = *b; // Set the first node to the second node
    *b = t; // Set the second node to the temporary variable
}

// A standard function to heapify at given idx. This function also updates position of nodes when they are swapped
void minHeapify(struct MinHeap *minHeap, int idx)
{
    int smallest, left, right; // Initialize the smallest, left, and right variables
    smallest = idx; // Initialize smallest as the root
    left = 2 * idx + 1; // Initialize left as the left child
    right = 2 * idx + 2; // Initialize right as the right child

    if (left < minHeap->size && minHeap->array[left]->dist < minHeap->array[smallest]->dist)
        smallest = left; // If the left child is smaller than the root, set smallest to the left child

    if (right < minHeap->size && minHeap->array[right]->dist < minHeap->array[smallest]->dist)
        smallest = right; // If the right child is smaller than the root, set smallest to the right child

    if (smallest != idx) // If the smallest node is not the root
    {
        // The nodes to be swapped in min heap
        struct MinHeapNode *smallestNode = minHeap->array[smallest]; // Store the smallest node
        struct MinHeapNode *idxNode = minHeap->array[idx]; // Store the root node

        // Swap positions
        minHeap->pos[smallestNode->v] = idx; // Set the position of the smallest node to the root
        minHeap->pos[idxNode->v] = smallest; // Set the position of the root to the smallest node

        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]); // Swap the smallest node and the root

        minHeapify(minHeap, smallest); // Heapify the smallest node
    }
}

// A utility function to check if the given minHeap is empty or not
int isEmpty(struct MinHeap *minHeap)
{
    return minHeap->size == 0; // Return 1 if the heap is empty, 0 otherwise
}

// Standard function to extract minimum node from heap
struct MinHeapNode *extractMin(struct MinHeap *minHeap)
{
    if (isEmpty(minHeap)) // If the heap is empty
        return NULL; // Return NULL

    // Store the root node
    struct MinHeapNode *root = minHeap->array[0]; // Store the root node

    // Replace root node with last node
    struct MinHeapNode *lastNode = minHeap->array[minHeap->size - 1]; // Store the last node
    minHeap->array[0] = lastNode; // Set the root node to the last node

    // Update position of last node
    minHeap->pos[root->v] = minHeap->size - 1; // Set the position of the root node to the last node
    minHeap->pos[lastNode->v] = 0; // Set the position of the last node to the root node

    // Reduce heap size and heapify root
    --minHeap->size; // Decrement the size of the heap
    minHeapify(minHeap, 0); // Heapify the root node

    return root;
}

// Function to decrease dist value of a given vertex v. This function uses pos[] of min heap to get the current index of node in min heap
void decreaseKey(struct MinHeap *minHeap, int v, int dist) {
    int i = minHeap->pos[v]; // Get the index of v in the heap array
    minHeap->array[i]->dist = dist; // Set the distance of v to dist

    // Optimized loop
    while (i && minHeap->array[i]->dist < minHeap->array[(i - 1) / 2]->dist) {
        minHeap->pos[minHeap->array[i]->v] = (i - 1) / 2; // Set the position of the node at i to the parent of the node at i
        minHeap->pos[minHeap->array[(i - 1) / 2]->v] = i; // Set the position of the parent of the node at i to i
        i = (i - 1) / 2; // Set i to the parent of i
    }
    swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]); // Swap the node at i and the parent of the node at i
}

// A utility function to print the path from source to j using parent array
void printPath(int parent[], int j)
{
    // Base Case : If j is source
    if (parent[j] == -1) // If the parent of j is -1
        return;

    printPath(parent, parent[j]); // Recursively print the path from source to parent[j]

    printf(" -> %d", j); // Print the vertex
}

// A utility function to print the constructed distance array
void printArr(int dist[], int n, int parent[], int src)
{
    printf("Vertex\t Distance\tPath"); // Print the column headers
    for (int i = 0; i < n; i++) // For each vertex
    {
        if (i != src && dist[i] != INT_MAX) // If the vertex is not the source and the distance is not infinity
        {
            printf("\n%d -> %d \t\t %d\t\t%d", src, i, dist[i], src); // Print the source, vertex, distance, and source
            printPath(parent, i); // Print the path from source to i
        }
    }
}

// A utility function to print the adjacency matrix
void printMatrix(int **graph, int vertices)
{
    printf("Adjacency Matrix:\n");

    // Print the column headers
    printf("    ");
    for (int j = 0; j < vertices; j++) // For each vertex
    {
        printf("%4d", j); // Print the vertex
    }
    printf("\n");

    // Print the row headers and the matrix values
    for (int i = 0; i < vertices; i++) // For each vertex
    {
        printf("%4d", i); // Print the vertex
        for (int j = 0; j < vertices; j++) // For each vertex
        {
            if (graph[i][j] == INT_MAX) // If the distance is infinity
            {
                printf("%4s", "INF"); // Print INF
            }
            else
            {
                printf("%4d", graph[i][j]); // Print the distance
            }
        }
        printf("\n");
    }
}

void dijkstra(int **graph, int src, int vertices)
{
    int dist[vertices]; // The output array. dist[i] will hold the shortest distance from src to i
    int parent[vertices]; // Parent array to store shortest path tree

    // minHeap represents set E
    struct MinHeap *minHeap = createMinHeap(vertices); // Create the min heap

    // Initialize min heap with all vertices. dist value of all vertices
    for (int v = 0; v < vertices; ++v) // For each vertex
    {
        parent[src] = -1; // Set the parent of the source to -1
        dist[v] = INT_MAX; // Set the distance to infinity
        minHeap->array[v] = newMinHeapNode(v, dist[v]); // Create a new min heap node
        minHeap->pos[v] = v; // Set the position of the node to v
    }

    // Make dist value of src vertex as 0 so that it is extracted first
    minHeap->array[src] = newMinHeapNode(src, dist[src]); // Create a new min heap node
    minHeap->pos[src] = src; // Set the position of the node to src
    dist[src] = 0; // Set the distance to 0
    decreaseKey(minHeap, src, dist[src]); // Decrease the key of the node at src

    // Initially size of min heap is equal to V
    minHeap->size = vertices; // Set the size of the heap to the number of vertices

    // In the followin loop, min heap contains all nodes whose shortest distance is not yet finalized.
    while (!isEmpty(minHeap)) // While the heap is not empty
    {
        // Extract the vertex with minimum distance value
        struct MinHeapNode *minHeapNode = extractMin(minHeap); // Extract the minimum node
        int u = minHeapNode->v; // Store the extracted vertex number

        // Traverse through all adjacent vertices of u (the extracted vertex) and update their distance values
        for (int v = 0; v < vertices; ++v) // For each vertex
        {
            // If shortest distance to v is not finalized yet, and distance to v through u is less than its previously calculated distance
            if (graph[u][v] != INT_MAX && dist[u] != INT_MAX && dist[u] + graph[u][v] < dist[v]) // If the distance from u to v is not infinity and the distance from the source to u is not infinity and the distance from the source to v through u is less than the distance from the source to v
            {
                dist[v] = dist[u] + graph[u][v]; // Set the distance from the source to v to the distance from the source to u plus the distance from u to v
                parent[v] = u; // Set the parent of v to u
                decreaseKey(minHeap, v, dist[v]); // Decrease the key of v
            }
        }
    }

    // Print the calculated shortest distances
    printArr(dist, vertices, parent, src);

    // Free the allocated memory
    free(minHeap->pos); // Free the position array
    free(minHeap->array); // Free the array
    free(minHeap); // Free the heap

    printf("\n");

    return;
}

// Function to read the graph from file
void readGraphFromFile(int **graph, const char *filename)
{
    FILE *file = fopen(filename, "r"); // Open the file for reading
    if (file == NULL) // If the file cannot be opened
    {
        perror("Error opening file"); // Print an error message
        return;
    }
    int src, dest, weight; // Initialize the source, destination, and weight variables
    while (fscanf(file, "%d %d %d", &src, &dest, &weight) != EOF) // While there is another line in the file
    {
        graph[src][dest] = weight; // Set the distance from the source to the destination to the weight
        graph[dest][src] = weight; // Set the distance from the destination to the source to the weight
    }
    fclose(file); // Close the file
}

// Function to get the number of vertices in the graph
int getNumberOfVertices(const char *filename) {
    FILE *file = fopen(filename, "r"); // Open the file for reading
    if (file == NULL) { // If the file cannot be opened
        perror("Error opening file"); // Print an error message
        return -1; // Return -1 to indicate an error
    }

    int src, dest, weight; // Initialize the source, destination, and weight variables
    int maxVertexIndex = -1; // Initialize the maximum vertex index to -1
    while (fscanf(file, "%d %d %d", &src, &dest, &weight) != EOF) { // While there is another line in the file
        if (src > maxVertexIndex) maxVertexIndex = src; // If the source is greater than the maximum vertex index, set the maximum vertex index to the source
        if (dest > maxVertexIndex) maxVertexIndex = dest; // If the destination is greater than the maximum vertex index, set the maximum vertex index to the destination
    }

    fclose(file); // Close the file

    // The number of vertices is 1 more than the maximum index since vertices start at 0
    return maxVertexIndex + 1; // Return the number of vertices
}

int main()
{
    const char *filename = "graph.txt"; // Initialize the filename variable
    int vertices = getNumberOfVertices(filename); // Get the number of vertices in the graph
    if (vertices == -1) { // If there was an error getting the number of vertices
        printf("Error reading graph file\n"); // Print an error message
        return 1; // Return 1 to indicate an error
    }

    int **graph = createGraph(vertices); // Create the graph

    // Read graph from file
    readGraphFromFile(graph, filename);

    // Print the adjacency matrix
    printMatrix(graph, vertices);

    // Run Dijkstra's algorithm with user input
    int src; // Initialize the source variable
    printf("Enter the source vertex: "); // Prompt the user for the source vertex
    scanf("%d", &src); // Read in the source vertex
    dijkstra(graph, src, vertices); // Run Dijkstra's algorithm

    // Free the allocated memory
    for (int i = 0; i < vertices; i++) // For each vertex
    {
        free(graph[i]); // Free the row
    }
    free(graph); // Free the graph

    return 0;
}