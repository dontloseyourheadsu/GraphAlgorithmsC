#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define MAX_SIZE 100

/*GRAPH*/

typedef struct edge {
    int vertex;         // Vertex index
    int weight;         // Weight of the edge
    struct edge *next;  // Pointer to the next edge
} Edge;

typedef struct {
    Edge *head;  // Head pointer for the list of edges
} Vertex;

typedef struct {
    int numVertices;  // Number of vertices in the graph
    Vertex *vertices; // Array of vertices
} Graph;

Graph* createGraph(int numVertices) {
    Graph *graph = (Graph *)malloc(sizeof(Graph));
    graph->numVertices = numVertices;
    graph->vertices = (Vertex *)malloc(numVertices * sizeof(Vertex));

    for (int i = 0; i < numVertices; i++) {
        graph->vertices[i].head = NULL;
    }

    return graph;
}

void addEdge(Graph *graph, int src, int dest, int weight) {
    Edge *newEdge = (Edge *)malloc(sizeof(Edge));
    newEdge->vertex = dest;
    newEdge->weight = weight;
    newEdge->next = graph->vertices[src].head;
    graph->vertices[src].head = newEdge;
}

void printGraph(Graph *graph) {
    for (int i = 0; i < graph->numVertices; i++) {
        Edge *curr = graph->vertices[i].head;
        printf("Vertex %d: ", i);
        while (curr != NULL) {
            printf("(%d, %d) ", curr->vertex, curr->weight);
            curr = curr->next;
        }
        printf("\n");
    }
}

int getNumVertices(Graph *graph) {
    return graph->numVertices;
}

/*PRIORITY QUEUE (MIN HEAP)*/

typedef struct {
    int vertex;  // Vertex index
    int key;     // Key value of the vertex
} HeapNode;

typedef struct {
    int capacity;    // Capacity of the heap
    int size;        // Size of the heap
    int *pos;        // Position of the vertex in the heap
    HeapNode **heap; // Array of heap nodes
} PriorityQueue;

PriorityQueue* createPriorityQueue(int capacity) {
    PriorityQueue *pq = (PriorityQueue *)malloc(sizeof(PriorityQueue));
    pq->capacity = capacity;
    pq->size = 0;
    pq->pos = (int *)malloc(capacity * sizeof(int));
    pq->heap = (HeapNode **)malloc(capacity * sizeof(HeapNode *));
    return pq;
}

void swapHeapNodes(HeapNode **a, HeapNode **b) {
    HeapNode *temp = *a;
    *a = *b;
    *b = temp;
}

void minHeapify(PriorityQueue *pq, int idx) {
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < pq->size && pq->heap[left]->key < pq->heap[smallest]->key) {
        smallest = left;
    }

    if (right < pq->size && pq->heap[right]->key < pq->heap[smallest]->key) {
        smallest = right;
    }

    if (smallest != idx) {
        HeapNode *smallestNode = pq->heap[smallest];
        HeapNode *idxNode = pq->heap[idx];

        pq->pos[smallestNode->vertex] = idx;
        pq->pos[idxNode->vertex] = smallest;

        swapHeapNodes(&pq->heap[smallest], &pq->heap[idx]);

        minHeapify(pq, smallest);
    }
}

int isEmpty(PriorityQueue *pq) {
    return pq->size == 0;
}

HeapNode* extractMin(PriorityQueue *pq) {
    if (isEmpty(pq)) {
        return NULL;
    }

    HeapNode *root = pq->heap[0];
    HeapNode *lastNode = pq->heap[pq->size - 1];
    pq->heap[0] = lastNode;

    pq->pos[root->vertex] = pq->size - 1;
    pq->pos[lastNode->vertex] = 0;

    pq->size--;
    minHeapify(pq, 0);

    return root;
}

void decreaseKey(PriorityQueue *pq, int vertex, int key) {
    int i = pq->pos[vertex];
    pq->heap[i]->key = key;

    while (i && pq->heap[i]->key < pq->heap[(i - 1) / 2]->key) {
        pq->pos[pq->heap[i]->vertex] = (i - 1) / 2;
        pq->pos[pq->heap[(i - 1) / 2]->vertex] = i;
        swapHeapNodes(&pq->heap[i], &pq->heap[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

int isInPriorityQueue(PriorityQueue *pq, int vertex) {
    return pq->pos[vertex] < pq->size;
}

/*DIJKSTRA ALGORITHM (DISTANCE AND PATH)*/
void printPath(int parent[], int j) {
    // Create an empty stack
    int stack[MAX_SIZE];
    int top = -1;

    // Push the destination node onto the stack
    stack[++top] = j;

    // While the top of the stack is not the source node
    while (parent[stack[top]] != -1) {
        // Push the parent of the top of the stack onto the stack
        stack[++top] = parent[stack[top]];
    }

    // Pop elements from the stack and print them until the stack is empty
    while (top != -1) {
        printf("%d ", stack[top--]);
    }
}

void printSolution(int dist[], int parent[], int src, int dest) {
    printf("Vertex\t Distance\tPath");
    printf("\n%d -> %d \t\t %d\t\t%d ", src, dest, dist[dest], src);
    printPath(parent, dest);
    printf("\n");
}

void Dijkstra(Graph *graph, int src) {
    int numVertices = graph->numVertices;
    int dist[numVertices];
    int parent[numVertices];

    PriorityQueue *pq = createPriorityQueue(numVertices);

    for (int v = 0; v < numVertices; v++) {
        dist[v] = INT_MAX;
        parent[v] = -1;
        pq->heap[v] = (HeapNode *)malloc(sizeof(HeapNode));
        pq->heap[v]->vertex = v;
        pq->heap[v]->key = dist[v];
        pq->pos[v] = v;
    }

    pq->heap[src] = (HeapNode *)malloc(sizeof(HeapNode));
    pq->heap[src]->vertex = src;
    pq->heap[src]->key = dist[src];
    pq->pos[src] = src;
    dist[src] = 0;
    decreaseKey(pq, src, dist[src]);

    pq->size = numVertices;

    while (!isEmpty(pq)) {
        HeapNode *minHeapNode = extractMin(pq);
        int u = minHeapNode->vertex;

        Edge *curr = graph->vertices[u].head;
        while (curr != NULL) {
            int v = curr->vertex;

            if (isInPriorityQueue(pq, v) && dist[u] != INT_MAX && curr->weight + dist[u] < dist[v]) {
                dist[v] = dist[u] + curr->weight;
                parent[v] = u;
                decreaseKey(pq, v, dist[v]);
            }
            curr = curr->next;
        }
    }

    printSolution(dist, parent, src, 6);
}

/*MAIN FUNCTION (FILE INPUT)*/
int main() {
    FILE *fp;
    fp = fopen("graph2.txt", "r");

    int numVertices;
    fscanf(fp, "%d", &numVertices);

    Graph *graph = createGraph(numVertices);

    int src, dest, weight;
    while (fscanf(fp, "%d %d %d", &src, &dest, &weight) != EOF) {
        addEdge(graph, src, dest, weight);
    }

    printGraph(graph);

    Dijkstra(graph, 0);

    return 0;
}
