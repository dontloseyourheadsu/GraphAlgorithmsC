#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define INITIAL_SIZE 5

typedef struct
{
    int vertex;
    int weight;
} Edge;

typedef struct
{
    Edge *edges;
    int edgeCount;
    int edgeSize;
} Node;

typedef struct
{
    Node *nodes;
    int nodeCount;
    int nodeSize;
} Graph;

Graph *createGraph(int initialSize)
{
    Graph *graph = (Graph *)malloc(sizeof(Graph));
    graph->nodes = (Node *)malloc(initialSize * sizeof(Node));
    graph->nodeCount = 0;
    graph->nodeSize = initialSize;
    for (int i = 0; i < initialSize; i++)
    {
        graph->nodes[i].edges = (Edge *)malloc(INITIAL_SIZE * sizeof(Edge));
        graph->nodes[i].edgeCount = 0;
        graph->nodes[i].edgeSize = INITIAL_SIZE;
    }
    return graph;
}

void addEdge(Graph *graph, int src, int dest, int weight)
{
    // Check if need to expand the array of edges
    if (graph->nodes[src].edgeCount >= graph->nodes[src].edgeSize)
    {
        graph->nodes[src].edgeSize *= 2;
        graph->nodes[src].edges = (Edge *)realloc(graph->nodes[src].edges, graph->nodes[src].edgeSize * sizeof(Edge));
    }
    graph->nodes[src].edges[graph->nodes[src].edgeCount].vertex = dest;
    graph->nodes[src].edges[graph->nodes[src].edgeCount].weight = weight;
    graph->nodes[src].edgeCount++;
}

void freeGraph(Graph *graph)
{
    for (int i = 0; i < graph->nodeSize; i++)
    {
        free(graph->nodes[i].edges);
    }
    free(graph->nodes);
    free(graph);
}

typedef struct
{
    int vertex;
    int distance;
} MinHeapNode;

typedef struct
{
    MinHeapNode *elements;
    int size;
    int capacity;
} MinHeap;

MinHeap *createMinHeap(int capacity)
{
    MinHeap *minHeap = (MinHeap *)malloc(sizeof(MinHeap));
    minHeap->elements = (MinHeapNode *)malloc(capacity * sizeof(MinHeapNode));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    return minHeap;
}

void swapMinHeapNode(MinHeapNode *a, MinHeapNode *b)
{
    MinHeapNode temp = *a;
    *a = *b;
    *b = temp;
}

void minHeapify(MinHeap *minHeap, int idx)
{
    int smallest = idx;
    int left = 2 * idx + 1;
    int right = 2 * idx + 2;

    if (left < minHeap->size && minHeap->elements[left].distance < minHeap->elements[smallest].distance)
    {
        smallest = left;
    }

    if (right < minHeap->size && minHeap->elements[right].distance < minHeap->elements[smallest].distance)
    {
        smallest = right;
    }

    if (smallest != idx)
    {
        swapMinHeapNode(&minHeap->elements[smallest], &minHeap->elements[idx]);
        minHeapify(minHeap, smallest);
    }
}

MinHeapNode extractMin(MinHeap *minHeap)
{
    if (minHeap->size == 0)
    {
        return (MinHeapNode){-1, 0};
    }

    MinHeapNode root = minHeap->elements[0];
    minHeap->elements[0] = minHeap->elements[minHeap->size - 1];
    minHeap->size--;
    minHeapify(minHeap, 0);

    return root;
}

void decreaseKey(MinHeap *minHeap, int vertex, int distance)
{
    int i;
    for (i = 0; i < minHeap->size; ++i)
    {
        if (minHeap->elements[i].vertex == vertex)
        {
            break;
        }
    }

    minHeap->elements[i].distance = distance;

    while (i != 0 && minHeap->elements[(i - 1) / 2].distance > minHeap->elements[i].distance)
    {
        swapMinHeapNode(&minHeap->elements[i], &minHeap->elements[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

void insertMinHeap(MinHeap *minHeap, MinHeapNode minHeapNode)
{
    if (minHeap->size == minHeap->capacity)
    {
        return;
    }

    minHeap->size++;
    int i = minHeap->size - 1;
    minHeap->elements[i] = minHeapNode;

    while (i != 0 && minHeap->elements[(i - 1) / 2].distance > minHeap->elements[i].distance)
    {
        swapMinHeapNode(&minHeap->elements[i], &minHeap->elements[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

void freeMinHeap(MinHeap *minHeap)
{
    free(minHeap->elements);
    free(minHeap);
}

void dijkstra(Graph *graph, int src)
{
    int V = graph->nodeCount;
    int dist[V];
    int prev[V];

    MinHeap *minHeap = createMinHeap(V);

    for (int v = 0; v < V; ++v)
    {
        dist[v] = INT_MAX;
        prev[v] = -1;
        minHeap->elements[v].vertex = v;
        minHeap->elements[v].distance = dist[v];
        minHeap->size++;
    }

    dist[src] = 0;
    decreaseKey(minHeap, src, dist[src]);

    while (minHeap->size != 0)
    {
        MinHeapNode minHeapNode = extractMin(minHeap);
        int u = minHeapNode.vertex;

        for (int i = 0; i < graph->nodes[u].edgeCount; ++i)
        {
            int v = graph->nodes[u].edges[i].vertex;

            if (isInMinHeap(minHeap, v) && dist[u] != INT_MAX && graph->nodes[u].edges[i].weight + dist[u] < dist[v])
            {
                dist[v] = dist[u] + graph->nodes[u].edges[i].weight;
                prev[v] = u;
                decreaseKey(minHeap, v, dist[v]);
            }
        }
    }

    printf("Vertex\tDistance\tPath\n");
    for (int i = 0; i < V; ++i)
    {
        printf("%d\t%d\t\t", i, dist[i]);
        int j = i;
        while (j != -1)
        {
            printf("%d ", j);
            j = prev[j];
        }
        printf("\n");
    }

    freeMinHeap(minHeap);
}

int isInMinHeap(MinHeap *minHeap, int v)
{
    for (int i = 0; i < minHeap->size; ++i)
    {
        if (minHeap->elements[i].vertex == v)
        {
            return 1;
        }
    }
    return 0;
}

void readGraphFromFile(char *filename, Graph *graph)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        perror("Error opening file");
        return;
    }

    int src, dest, weight;
    int maxNodeIndex = -1; // Keep track of the highest node index

    while (fscanf(file, "%d %d %d", &src, &dest, &weight) != EOF)
    {
        if (src >= graph->nodeSize || dest >= graph->nodeSize)
        {
            // Logic to expand the graph size (if needed)...
        }
        addEdge(graph, src, dest, weight);
        if (src > maxNodeIndex)
            maxNodeIndex = src;
        if (dest > maxNodeIndex)
            maxNodeIndex = dest;
    }

    graph->nodeCount = maxNodeIndex + 1; // Update node count

    fclose(file);
}

void printGraph(Graph *graph)
{
    int V = graph->nodeSize;
    int **adjMatrix = (int **)malloc(V * sizeof(int *));
    for (int i = 0; i < V; i++)
    {
        adjMatrix[i] = (int *)malloc(V * sizeof(int));
        for (int j = 0; j < V; j++)
        {
            adjMatrix[i][j] = 0; // Initialize with 0
        }
    }

    // Fill the adjacency matrix
    for (int i = 0; i < V; i++)
    {
        for (int j = 0; j < graph->nodes[i].edgeCount; j++)
        {
            int dest = graph->nodes[i].edges[j].vertex;
            adjMatrix[i][dest] = graph->nodes[i].edges[j].weight;
        }
    }

    // Print the adjacency matrix
    printf("Adjacency Matrix:\n");

    // Print column labels
    printf("   ");
    for (int j = 0; j < V; j++)
    {
        printf("%d ", j);
    }
    printf("\n");

    for (int i = 0; i < V; i++)
    {
        // Print row label
        printf("%d: ", i);
        for (int j = 0; j < V; j++)
        {
            printf("%d ", adjMatrix[i][j]);
        }
        printf("\n");
    }

    // Free the allocated memory
    for (int i = 0; i < V; i++)
    {
        free(adjMatrix[i]);
    }
    free(adjMatrix);
}

int main()
{
    char filename[100] = "graph.txt";

    Graph *graph = createGraph(INITIAL_SIZE);
    readGraphFromFile(filename, graph);

    printGraph(graph);

    int modify;
    printf("Do you want to modify the graph? (1 for Yes, 0 for No): ");
    scanf("%d", &modify);

    if (modify)
    {
        int src, dest, weight;
        printf("Enter source, destination, and weight for new edge: ");
        scanf("%d %d %d", &src, &dest, &weight);
        addEdge(graph, src, dest, weight);
        printGraph(graph);
    }

    int source;
    printf("Enter the source vertex for Dijkstra's algorithm: ");
    scanf("%d", &source);

    dijkstra(graph, source);

    freeGraph(graph);
    return 0;
}
