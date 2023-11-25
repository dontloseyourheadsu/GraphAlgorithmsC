#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define MAX_VERTICES 100

struct MinHeapNode
{
    int v;
    int dist;
};

struct MinHeap
{
    int size;
    int capacity;
    struct MinHeapNode **array;
    int *pos; // This is needed for decreaseKey()
};

struct MinHeapNode *newMinHeapNode(int v, int dist)
{
    struct MinHeapNode *minHeapNode = (struct MinHeapNode *)malloc(sizeof(struct MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->dist = dist;
    return minHeapNode;
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

struct MinHeap *createMinHeap(int capacity)
{
    struct MinHeap *minHeap = (struct MinHeap *)malloc(sizeof(struct MinHeap));
    minHeap->pos = (int *)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = (struct MinHeapNode **)malloc(capacity * sizeof(struct MinHeapNode *));
    return minHeap;
}

void swapMinHeapNode(struct MinHeapNode **a, struct MinHeapNode **b)
{
    struct MinHeapNode *t = *a;
    *a = *b;
    *b = t;
}

void minHeapify(struct MinHeap *minHeap, int idx)
{
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < minHeap->size && minHeap->array[left]->dist < minHeap->array[smallest]->dist)
        smallest = left;

    if (right < minHeap->size && minHeap->array[right]->dist < minHeap->array[smallest]->dist)
        smallest = right;

    if (smallest != idx)
    {
        // The nodes to be swapped in min heap
        struct MinHeapNode *smallestNode = minHeap->array[smallest];
        struct MinHeapNode *idxNode = minHeap->array[idx];

        // Swap positions
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;

        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);

        minHeapify(minHeap, smallest);
    }
}

int isEmpty(struct MinHeap *minHeap)
{
    return minHeap->size == 0;
}

struct MinHeapNode *extractMin(struct MinHeap *minHeap)
{
    if (isEmpty(minHeap))
        return NULL;

    // Store the root node
    struct MinHeapNode *root = minHeap->array[0];

    // Replace root node with last node
    struct MinHeapNode *lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    // Update position of last node
    minHeap->pos[root->v] = minHeap->size - 1;
    minHeap->pos[lastNode->v] = 0;

    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);

    return root;
}

void decreaseKey(struct MinHeap *minHeap, int v, int dist) {
    int i = minHeap->pos[v];
    minHeap->array[i]->dist = dist;

    // Optimized loop
    while (i && minHeap->array[i]->dist < minHeap->array[(i - 1) / 2]->dist) {
        minHeap->pos[minHeap->array[i]->v] = (i - 1) / 2;
        minHeap->pos[minHeap->array[(i - 1) / 2]->v] = i;
        i = (i - 1) / 2;
    }
    swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);
}

// A utility function to print the path from source to j using parent array
void printPath(int parent[], int j)
{
    // Base Case : If j is source
    if (parent[j] == -1)
        return;

    printPath(parent, parent[j]);

    printf(" -> %d", j);
}

// A utility function to print the constructed distance array
void printArr(int dist[], int n, int parent[], int src)
{
    printf("Vertex\t Distance\tPath");
    for (int i = 0; i < n; i++)
    {
        if (i != src && dist[i] != INT_MAX)
        {
            printf("\n%d -> %d \t\t %d\t\t%d", src, i, dist[i], src);
            printPath(parent, i); // Print the path from source to i
        }
    }
}

void printMatrix(int **graph, int vertices)
{
    printf("Adjacency Matrix:\n");

    // Print the column headers
    printf("    ");
    for (int j = 0; j < vertices; j++)
    {
        printf("%4d", j);
    }
    printf("\n");

    // Print the row headers and the matrix values
    for (int i = 0; i < vertices; i++)
    {
        printf("%4d", i);
        for (int j = 0; j < vertices; j++)
        {
            if (graph[i][j] == INT_MAX)
            {
                printf("%4s", "INF");
            }
            else
            {
                printf("%4d", graph[i][j]);
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
    struct MinHeap *minHeap = createMinHeap(vertices);

    // Initialize min heap with all vertices. dist value of all vertices
    for (int v = 0; v < vertices; ++v)
    {
        parent[src] = -1;
        dist[v] = INT_MAX;
        minHeap->array[v] = newMinHeapNode(v, dist[v]);
        minHeap->pos[v] = v;
    }

    // Make dist value of src vertex as 0 so that it is extracted first
    minHeap->array[src] = newMinHeapNode(src, dist[src]);
    minHeap->pos[src] = src;
    dist[src] = 0;
    decreaseKey(minHeap, src, dist[src]);

    // Initially size of min heap is equal to V
    minHeap->size = vertices;

    // In the followin loop, min heap contains all nodes whose shortest distance is not yet finalized.
    while (!isEmpty(minHeap))
    {
        // Extract the vertex with minimum distance value
        struct MinHeapNode *minHeapNode = extractMin(minHeap);
        int u = minHeapNode->v; // Store the extracted vertex number

        // Traverse through all adjacent vertices of u (the extracted vertex) and update their distance values
        for (int v = 0; v < vertices; ++v)
        {
            // If shortest distance to v is not finalized yet, and distance to v through u is less than its previously calculated distance
            if (graph[u][v] != INT_MAX && dist[u] != INT_MAX && dist[u] + graph[u][v] < dist[v])
            {
                dist[v] = dist[u] + graph[u][v];
                parent[v] = u;
                decreaseKey(minHeap, v, dist[v]);
            }
        }
    }

    // Print the calculated shortest distances
    printArr(dist, vertices, parent, src);

    // Free the allocated memory
    free(minHeap->pos);
    free(minHeap->array);
    free(minHeap);

    printf("\n");

    return;
}

// Function to read the graph from file
void readGraphFromFile(int **graph, const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        perror("Error opening file");
        return;
    }
    int src, dest, weight;
    while (fscanf(file, "%d %d %d", &src, &dest, &weight) != EOF)
    {
        graph[src][dest] = weight;
        graph[dest][src] = weight; // If the graph is undirected
    }
    fclose(file);
}

int getNumberOfVertices(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return -1; // Return -1 to indicate an error
    }

    int src, dest, weight;
    int maxVertexIndex = -1;
    while (fscanf(file, "%d %d %d", &src, &dest, &weight) != EOF) {
        if (src > maxVertexIndex) maxVertexIndex = src;
        if (dest > maxVertexIndex) maxVertexIndex = dest;
    }

    fclose(file);

    // The number of vertices is 1 more than the maximum index since vertices start at 0
    return maxVertexIndex + 1;
}

int main()
{
    // Determine the number of vertices by reading the graph file
    const char *filename = "graph.txt";
    int vertices = getNumberOfVertices(filename);
    if (vertices == -1) {
        printf("Error reading graph file\n");
        return 1;
    }

    int **graph = createGraph(vertices);

    // Read graph from file
    readGraphFromFile(graph, filename);

    // Print the adjacency matrix
    printMatrix(graph, vertices);

    // Run Dijkstra's algorithm with user input
    int src;
    printf("Enter the source vertex: ");
    scanf("%d", &src);
    dijkstra(graph, src, vertices);

    // Free the allocated memory
    for (int i = 0; i < vertices; i++)
    {
        free(graph[i]);
    }
    free(graph);

    return 0;
}