#include "btree_mgr.h"
#include "buffer_mgr.h"
#include "storage_mgr.h"
#include "dberror.h"
#include "tables.h"
#include <string.h>
#include <stdlib.h>

typedef struct BTreeData {
	BM_BufferPool* bp;
	int n;
	int treeHeight;
	int lastPage;
	int entries;
	int rootPage;
	int* walkTreeIndex;
}BTreeData;

typedef struct BTreeNode {
	bool leaf;
	int selfPage;
	int numKeys;
	int parentPage;
	int leftPage;
	int rightPage;
	char** keyData;
	int* children;
	RID* data;
} BTreeNode;

typedef struct TreeScanData {
	int pageNum;
	int index;
	BTreeNode* BTN;
} TreeScanData;

RC readBTreeNode(BTreeHandle* tree, int pageNum, BTreeNode** BTN) {
	BTreeData* BTD = (BTreeData*)tree->mgmtData;
	BM_PageHandle* page = malloc(sizeof(BM_PageHandle));

	if (!page) {
		return RC_WRITE_FAILED;
	}

	RC rc = pinPage(BTD->bp, page, pageNum);
	if (rc != RC_OK) {
		free(page);
		return rc;
	}

	(*BTN) = malloc(sizeof(BTreeNode));
	if (!(*BTN)) {
		unpinPage(BTD->bp, page);
		free(page);
		return RC_WRITE_FAILED; // Memory allocation failure
	}

	int offset = 0;
	memcpy(&((*BTN)->leaf), page->data + offset, sizeof(bool));
	offset += sizeof(bool);
	memcpy(&((*BTN)->numKeys), page->data + offset, sizeof(int));
	offset += sizeof(int);
	memcpy(&((*BTN)->parentPage), page->data + offset, sizeof(int));
	offset += sizeof(int);
	memcpy(&((*BTN)->leftPage), page->data + offset, sizeof(int));
	offset += sizeof(int);
	memcpy(&((*BTN)->rightPage), page->data + offset, sizeof(int));
	offset += sizeof(int);

	// Allocate memory for keyData as an array of char pointers (string keys)
	(*BTN)->keyData = (char**)malloc(sizeof(char*) * BTD->n);
	if ((*BTN)->keyData == NULL) {
		// Handle memory allocation failure
		free(*BTN);
		unpinPage(BTD->bp, page);
		free(page);
		return RC_WRITE_FAILED;
	}

	for (int i = 0; i < (*BTN)->numKeys; i++) {
		// Assuming a fixed size for each key, e.g., 20 characters
		(*BTN)->keyData[i] = (char*)malloc(sizeof(char) * 20);
		memcpy((*BTN)->keyData[i], page->data + offset, 20);
		offset += 20;
	}

	if ((*BTN)->leaf == TRUE) {
		(*BTN)->data = malloc(sizeof(RID) * BTD->n);
		for (int i = 0; i < (*BTN)->numKeys; i++) {
			memcpy(&(((*BTN)->data[i]).page), page->data + offset, sizeof(int));
			offset += sizeof(int);
			memcpy(&(((*BTN)->data[i]).slot), page->data + offset, sizeof(int));
			offset += sizeof(int);
		}
	}
	else {
		(*BTN)->children = malloc(sizeof(int) * (BTD->n + 1));
		for (int i = 0; i < (*BTN)->numKeys + 1; i++) {
			memcpy(&((*BTN)->children[i]), page->data + offset, sizeof(int));
			offset += sizeof(int);
		}
	}

	unpinPage(BTD->bp, page);
	free(page);
	return RC_OK;
}

RC writeBTreeNode(BTreeHandle* tree, int pageNum, BTreeNode* BTN) {
	BTreeData* BTD = (BTreeData*)tree->mgmtData;
	BM_PageHandle* page = malloc(sizeof(BM_PageHandle));
	if (!page)
		return RC_WRITE_FAILED;

	RC rc = pinPage(BTD->bp, page, pageNum);
	if (rc != RC_OK) {
		free(page);
		return rc;
	}

	markDirty(BTD->bp, page); // Page is going to be overwritten
	int offset = 0;
	memcpy(page->data + offset, &(BTN->leaf), sizeof(bool)); // leafType
	offset += sizeof(bool);
	memcpy(page->data + offset, &(BTN->numKeys), sizeof(int)); // numKeys
	offset += sizeof(int);
	memcpy(page->data + offset, &(BTN->parentPage), sizeof(int)); // parentPageNum
	offset += sizeof(int);
	memcpy(page->data + offset, &(BTN->leftPage), sizeof(int)); // leftSiblingPageNum
	offset += sizeof(int);
	memcpy(page->data + offset, &(BTN->rightPage), sizeof(int)); // rightSiblingPageNum
	offset += sizeof(int);

	for (int i = 0; i < BTN->numKeys; i++) {
		// Assuming a fixed size for each key, e.g., 20 characters
		memcpy(page->data + offset, BTN->keyData[i], 20); // string keys
		offset += 20;
	}

	if (BTN->leaf == TRUE) {
		for (int i = 0; i < BTN->numKeys; i++) {
			memcpy(page->data + offset, &((BTN->data[i]).page), sizeof(int)); // page num
			offset += sizeof(int);
			memcpy(page->data + offset, &((BTN->data[i]).slot), sizeof(int)); // slot num
			offset += sizeof(int);
		}
	}
	else {
		for (int i = 0; i < BTN->numKeys + 1; i++) {
			memcpy(page->data + offset, &(BTN->children[i]), sizeof(int)); // child page number
			offset += sizeof(int);
		}
	}

	rc = markDirty(BTD->bp, page);
	if (rc != RC_OK) {
		unpinPage(BTD->bp, page);
		free(page);
		return rc;
	}

	rc = unpinPage(BTD->bp, page);
	free(page);
	return rc;
}

RC createBTreeNode(BTreeHandle* tree, BTreeNode** BTN, bool isLeaf) {
	BTreeData* BTD = (BTreeData*)tree->mgmtData;

	// Allocate memory for the new B-tree node
	*BTN = (BTreeNode*)malloc(sizeof(BTreeNode));
	if (*BTN == NULL) {
		return RC_WRITE_FAILED; // Error code for memory allocation failure
	}

	// Initialize the B-tree node structure
	(*BTN)->leaf = isLeaf;
	(*BTN)->numKeys = 0;
	(*BTN)->parentPage = -1;
	(*BTN)->leftPage = -1;
	(*BTN)->rightPage = -1;
	(*BTN)->keyData = NULL;
	(*BTN)->children = NULL;
	(*BTN)->data = NULL;

	// Allocate memory for keyData as array of char pointers (string keys)
	(*BTN)->keyData = (char**)malloc(sizeof(char*) * BTD->n);
	if ((*BTN)->keyData == NULL) {
		free(*BTN);
		return RC_WRITE_FAILED; // Error code for memory allocation failure
	}
	// Initialize all key pointers to NULL
	for (int i = 0; i < BTD->n; ++i) {
		(*BTN)->keyData[i] = NULL;
	}

	if (isLeaf) {
		// Allocate memory for data if the node is a leaf
		(*BTN)->data = (RID*)malloc(sizeof(RID) * BTD->n);
		if ((*BTN)->data == NULL) {
			free((*BTN)->keyData);
			free(*BTN);
			return RC_WRITE_FAILED; // Error code for memory allocation failure
		}
	}
	else {
		// Allocate memory for children if the node is not a leaf
		(*BTN)->children = (int*)malloc(sizeof(int) * (BTD->n + 1));
		if ((*BTN)->children == NULL) {
			free((*BTN)->keyData);
			free(*BTN);
			return RC_WRITE_FAILED; // Error code for memory allocation failure
		}
	}

	return RC_OK;
}

RC freeBTreeNode(BTreeNode* BTN) {
	if (!BTN)
		return RC_OK;

	if (BTN->leaf && BTN->data != NULL) {
		free(BTN->data);
	}

	if (BTN->children != NULL) {
		free(BTN->children);
	}

	if (BTN->keyData != NULL) {
		free(BTN->keyData);
	}

	free(BTN);
	return RC_OK;
}

void printTreeNode(BTreeHandle* tree, int pageNum)
{
	BTreeData* BTD = (BTreeData*)tree->mgmtData;
	BTreeNode* BTN = 0;
	readBTreeNode(tree, pageNum, &BTN);
	
	printf("(%i)[", BTD->walkTreeIndex[pageNum]);
	if (BTN->leaf == TRUE)
	{
		for (int i = 0; i < BTN->numKeys; i++)
		{
			if (i > 0)
			{
				printf(",");
			}
			printf("%i.%i,%i", BTN->data[i].page, BTN->data[i].slot, BTN->keyData[i]);
		}
		if (BTN->rightPage == -1)
		{
			printf("]\n");
		}
		else
		{
			printf(",%i]\n", BTD->walkTreeIndex[BTN->rightPage]);
		}
	}
	else
	{
		int numChildNodes = BTN->numKeys + 1;
		int* childNodes = malloc(sizeof(int) * numChildNodes);
		for (int i = 0; i < BTN->numKeys; i++)
		{
			printf("%i,%i,", BTD->walkTreeIndex[BTN->children[i]], BTN->keyData[i]);
			childNodes[i] = BTN->children[i];
		}
		printf("%i]\n", BTD->walkTreeIndex[BTN->children[BTN->numKeys]]);
		childNodes[BTN->numKeys] = BTN->children[BTN->numKeys];
		freeBTreeNode(BTN);
		for (int i = 0; i < numChildNodes; i++)
		{
			printTreeNode(tree, childNodes[i]);
		}
	}
}
void buildTreeIndex(BTreeHandle* tree, int pageNum, int height, int* index)
{
	BTreeData* BTD = (BTreeData*)tree->mgmtData;
	BTreeNode* BTN = 0;
	int numChildNodes = 0;

	if (height+1 < BTD->treeHeight)
	{
		int* childNodes = malloc(sizeof(int) * (BTD->n + 1));
		readBTreeNode(tree, pageNum, &BTN);
		for (int i = 0; i < BTN->numKeys + 1; i++)
		{
			childNodes[i] = BTN->children[i];
		}
		numChildNodes = BTN->numKeys + 1;
		freeBTreeNode(BTN);
		for (int i = 0; i < numChildNodes; i++)
		{
			BTD->walkTreeIndex[childNodes[i]] = index[0];
			index[0]++;
			if (height + 2 < BTD->treeHeight)
			{
				buildTreeIndex(tree, childNodes[i], height + 1, index);
			}
		}
		free(childNodes);
	}
}

// init and shutdown index manager
extern RC initIndexManager(void* mgmtData)
{
	return RC_OK;
}
extern RC shutdownIndexManager()
{
	return RC_OK;
}

// create, destroy, open, and close an btree index
extern RC createBtree(char* idxId, DataType keyType, int n)
{
	SM_FileHandle fh;
	SM_PageHandle ph = (SM_PageHandle)malloc(PAGE_SIZE);
	createPageFile(idxId);
	openPageFile(idxId, &fh);
	readFirstBlock(&fh, ph);
	int treeHeight = 0;
	int offset = 0;
	int lastPage = 0;
	int rootPage = -1;
	int entries = 0;
	memcpy(ph + offset, &keyType, sizeof(DataType));//keyType
	offset += sizeof(DataType);
	memcpy(ph + offset, &n, sizeof(int));//n value for b+ tree
	offset += sizeof(int);
	memcpy(ph + offset, &rootPage, sizeof(int));//root Page number
	offset += sizeof(int);
	memcpy(ph + offset, &treeHeight, sizeof(int));//current height of tree
	offset += sizeof(int);
	memcpy(ph + offset, &lastPage, sizeof(int));//last used Page
	offset += sizeof(int);
	memcpy(ph + offset, &entries, sizeof(int));//number of entries
	offset += sizeof(int);
	writeCurrentBlock(&fh, ph);
	closePageFile(&fh);
	free(ph);
	return RC_OK;
}
extern RC openBtree(BTreeHandle** tree, char* idxId)
{
	*tree = 0;
	BM_BufferPool* bp = malloc(sizeof(BM_BufferPool));
	RC rc = initBufferPool(bp, idxId, 10, RS_FIFO, NULL);
	if (rc)// if file not found
	{
		free(bp);
		return rc;
	}
	*tree = malloc(sizeof(BTreeHandle));
	(*tree)->idxId = malloc(strlen(idxId)+1);
	strcpy((*tree)->idxId, idxId);

	BTreeData* BTD = malloc(sizeof(BTreeData));
	BTD->bp = bp;
	BM_PageHandle* page = malloc(sizeof(BM_PageHandle));
	pinPage(BTD->bp, page, 0);
	
	int offset = 0;
	memcpy(&((*tree)->keyType), page->data + offset, sizeof(DataType));//keyType
	offset += sizeof(DataType);
	memcpy(&(BTD->n), page->data + offset, sizeof(int));//n value for b+ tree
	offset += sizeof(int);
	memcpy(&(BTD->rootPage), page->data + offset, sizeof(int));//root Page number
	offset += sizeof(int);
	memcpy(&(BTD->treeHeight), page->data + offset, sizeof(int));//current height of tree
	offset += sizeof(int);
	memcpy(&(BTD->lastPage), page->data + offset, sizeof(int));//last used Page
	offset += sizeof(int);
	memcpy(&(BTD->entries), page->data + offset, sizeof(int));//number of entries
	offset += sizeof(int);
	unpinPage(BTD->bp, page);
	if (BTD->n > 5) //if tree is wider than 5, re-open buffer pool with larger pool.
	{
		shutdownBufferPool(bp);
		initBufferPool(bp, idxId, BTD->n * 2, RS_FIFO, NULL);
	}
	(*tree)->mgmtData = (void*)BTD;

	free(page);
	return RC_OK;
}
extern RC closeBtree(BTreeHandle* tree)
{
	if (tree == 0)
	{
		return RC_FILE_NOT_FOUND;
	}
	BTreeData* BTD = (BTreeData*)(tree->mgmtData);
	BM_PageHandle* page = malloc(sizeof(BM_PageHandle));
	RC rc = pinPage(BTD->bp, page, 0);
	if (rc)
	{
		free(page);
		return rc;
	}
	int offset = 0;
	memcpy(page->data + offset, &(tree->keyType), sizeof(DataType));//keyType
	offset += sizeof(DataType);
	memcpy(page->data + offset, &(BTD->n), sizeof(int));//n value for b+ tree
	offset += sizeof(int);
	memcpy(page->data + offset, &(BTD->rootPage), sizeof(int));//root Page number
	offset += sizeof(int);
	memcpy(page->data + offset, &(BTD->treeHeight), sizeof(int));//current height of tree
	offset += sizeof(int);
	memcpy( page->data + offset, &(BTD->lastPage), sizeof(int));//last used Page
	offset += sizeof(int);
	memcpy(page->data + offset, &(BTD->entries), sizeof(int));//number of entries
	offset += sizeof(int);
	markDirty(BTD->bp, page);
	unpinPage(BTD->bp, page);
	free(page);

	rc = forceFlushPool(BTD->bp);
	if (rc)
	{
		return rc;
	}
	shutdownBufferPool(BTD->bp);
	free(BTD->bp);
	free(BTD);
	free(tree->idxId);
	free(tree);
	tree = 0;
	return RC_OK;
}
extern RC deleteBtree(char* idxId)
{
	return destroyPageFile(idxId);
}

// access information about a b-tree
extern RC getNumNodes(BTreeHandle* tree, int* result)
{
	if (tree == 0)
	{
		return RC_FILE_NOT_FOUND;
	}
	BTreeData* BPD = (BTreeData*)tree->mgmtData;

	*result = BPD->lastPage;
	return RC_OK;
}
extern RC getNumEntries(BTreeHandle* tree, int* result)
{
	if (tree == 0)
	{
		return RC_FILE_NOT_FOUND;
	}
	BTreeData* BPD = (BTreeData*)tree->mgmtData;
	*result = BPD->entries;
	return RC_OK;
}
extern RC getKeyType(BTreeHandle* tree, DataType* result)
{
	if (tree == 0)
	{
		return RC_FILE_NOT_FOUND;
	}
	*result = tree->keyType;
	return RC_OK;
}

// index access
extern RC findKey(BTreeHandle* tree, Value* key, RID* result) {
	if (tree == NULL || key == NULL || result == NULL) {
		return RC_FILE_NOT_FOUND;
	}

	char* serializedSearchKey = serializeValue(key);  // Serialize the search key
	BTreeData* btreeData = (BTreeData*)tree->mgmtData;

	int currentPageNum = btreeData->rootPage;
	BTreeNode* currentNode = NULL;

	while (TRUE) {
		RC readStatus = readBTreeNode(tree, currentPageNum, &currentNode);
		if (readStatus != RC_OK) {
			free(serializedSearchKey);  // Free the serialized key on error
			return readStatus;
		}

		if (currentNode->leaf) {
			for (int i = 0; i < currentNode->numKeys; i++) {
				if (strcmp(currentNode->keyData[i], serializedSearchKey) == 0) {  // Compare serialized keys
					*result = currentNode->data[i];
					free(serializedSearchKey);
					free(currentNode);
					return RC_OK;
				}
			}
			free(serializedSearchKey);
			free(currentNode);
			return RC_IM_KEY_NOT_FOUND;  // Key not found in the leaf
		}
		else {
			// Traverse to find the appropriate child
			int i;
			for (i = 0; i < currentNode->numKeys; i++) {
				if (strcmp(serializedSearchKey, currentNode->keyData[i]) < 0) {  // Use serialized key for comparison
					break;
				}
			}
			currentPageNum = currentNode->children[i];
			free(currentNode);
		}
	}
}

extern RC insertKey(BTreeHandle* tree, Value* key, RID rid) {
	if (tree == NULL || key == NULL) {
		return RC_IM_KEY_NOT_FOUND;
	}

	char* serializedKey = serializeValue(key);
	BTreeData* treeData = (BTreeData*)tree->mgmtData;

	if (treeData->rootPage == -1) {
		// Create the first node as it's empty
		BTreeNode* newNode;
		RC createStatus = createBTreeNode(tree, &newNode, true); // true for leaf node
		if (createStatus != RC_OK) {
			free(serializedKey); // Free the serialized key
			return createStatus;
		}

		// Set the first available page number as the root page number
		newNode->selfPage = treeData->lastPage + 1;

		// Insert the serialized key into the new node
		newNode->keyData[0] = strdup(serializedKey); // Allocate and copy the serialized key
		newNode->data[0] = rid; // Set the corresponding RID
		newNode->numKeys = 1; // Update the number of keys in the new node

		// Update tree data
		treeData->rootPage = newNode->selfPage; // Use the assigned page number
		treeData->treeHeight = 1;
		treeData->entries = 1;

		// Write the new node to the disk
		RC writeStatus = writeBTreeNode(tree, newNode->selfPage, newNode);
		if (writeStatus != RC_OK) {
			freeBTreeNode(newNode); // Free the entire node structure if write fails
			free(serializedKey); // Free the serialized key
			return writeStatus;
		}

		treeData->lastPage = newNode->selfPage;
		freeBTreeNode(newNode); // Free the entire node structure when done
		free(serializedKey); // Free the serialized key
		return RC_OK;
	}

	BTreeNode* leafNode = NULL;

	// leafNode is found and assigned correctly
	RC findStatus = findLeafNode(tree, serializedKey, &leafNode);
	if (findStatus != RC_OK) {
		free(serializedKey);
		return findStatus;
	}

	if (leafNode == NULL) {
		free(serializedKey);
		return RC_IM_KEY_NOT_FOUND;
	}

	// Check for duplicates before trying to insert
	for (int i = 0; i < leafNode->numKeys; ++i) {
		if (strcmp(leafNode->keyData[i], serializedKey) == 0) {
			// Key already exists, update RID and avoid inserting a duplicate
			leafNode->data[i] = rid;
			free(serializedKey);
			RC writeStatus = writeBTreeNode(tree, leafNode->selfPage, leafNode);
			freeBTreeNode(leafNode); // Free the node
			return writeStatus;
		}
	}

	// Proceed with insertion if no duplicate was found
	if (leafNode->numKeys < treeData->n) {
		RC status = insertIntoLeaf(tree, leafNode, serializedKey, rid);
		free(serializedKey);
		freeBTreeNode(leafNode);
		return status;
	}
	else {
		char* newSerializedKey; // Placeholder for new key after split
		int newPageNum;
		RC status = splitLeafNode(tree, leafNode, serializedKey, rid, &newSerializedKey, &newPageNum);
		if (status != RC_OK) {
			free(serializedKey);
			freeBTreeNode(leafNode);
			return status;
		}
		status = insertIntoParent(tree, leafNode, newSerializedKey, newPageNum);
		free(serializedKey);
		free(newSerializedKey);
		freeBTreeNode(leafNode);
		return status;
	}
}
// helper-method for insertKey and deleteKey to find the leaf node for a given key
RC findLeafNode(BTreeHandle* tree, char* serializedKey, BTreeNode** leafNode) {
	// validating the input
	if (tree == NULL || leafNode == NULL || serializedKey == NULL) {
		return RC_IM_KEY_NOT_FOUND;
	}

	BTreeData* treeData = (BTreeData*)tree->mgmtData;
	int currentPageNum = treeData->rootPage; // Start from the root page
	BTreeNode* currentNode = NULL;

	// traversing from root to leaf
	while (TRUE) {
		// read the current node from page
		RC readStatus = readBTreeNode(tree, currentPageNum, &currentNode);
		if (readStatus != RC_OK) {
			// if read error
			if (currentNode != NULL) {
				freeBTreeNode(currentNode);
			}
			return readStatus;
		}

		if (currentNode->leaf) {
			*leafNode = currentNode;
			return RC_OK;
		}
		else {
			// traversing to the next node
			int i = 0;
			for (i = 0; i < currentNode->numKeys; i++) {
				if (strcmp(serializedKey, currentNode->keyData[i]) < 0) {
					break;
				}
			}
			// child node at index i
			currentPageNum = currentNode->children[i];

			// free the current node
			freeBTreeNode(currentNode);
		}
	}
	return RC_IM_KEY_NOT_FOUND;
}

// Helper-method for insertKey to insert a key-value pair into a leaf node
RC insertIntoLeaf(BTreeHandle* tree, BTreeNode* leafNode, char* serializedKey, RID rid) {
	// Validating the input parameters
	if (tree == NULL || leafNode == NULL || serializedKey == NULL) {
		return RC_IM_KEY_NOT_FOUND;
	}

	// Checking if leaf node is a leaf
	if (!leafNode->leaf) {
		return RC_IM_KEY_NOT_FOUND;
	}

	BTreeData* treeData = (BTreeData*)tree->mgmtData;

	// Checking if the leaf node has space to insert the key
	if (leafNode->numKeys < treeData->n) {
		// Inserting the key in the correct position to maintain order
		int i;
		for (i = leafNode->numKeys; i > 0 && strcmp(leafNode->keyData[i - 1], serializedKey) > 0; i--) {
			leafNode->keyData[i] = leafNode->keyData[i - 1]; // Shifting keys
			leafNode->data[i] = leafNode->data[i - 1]; // Shifting corresponding RIDs
		}

		// Allocate memory and copy the serialized key into the node
		leafNode->keyData[i] = strdup(serializedKey);
		leafNode->data[i] = rid;
		leafNode->numKeys++;

		// Writing the updated leaf node back to the disk
		return writeBTreeNode(tree, leafNode->selfPage, leafNode);
	}
	else {
		// The leaf node is full and needs to be split
		char* newSerializedKey; // Placeholder for the new key after splitting
		int newPageNum;

		// Split the full leaf node
		RC status = splitLeafNode(tree, leafNode, serializedKey, rid, &newSerializedKey, &newPageNum);
		if (status != RC_OK) {
			return status;
		}

		// The splitLeafNode function should return the new key and the new page number
		// Determine which node to insert into (the original or the new node after split)
		if (strcmp(serializedKey, newSerializedKey) < 0) {
			// Insert the key into the original node since the key is smaller
			return insertIntoLeaf(tree, leafNode, serializedKey, rid);
		}
		else {
			// Read the new leaf node from the disk
			BTreeNode* newLeafNode;
			status = readBTreeNode(tree, newPageNum, &newLeafNode);
			if (status != RC_OK) {
				return status;
			}

			// Insert the key into the new leaf node
			status = insertIntoLeaf(tree, newLeafNode, serializedKey, rid);

			// Write the new leaf node back to disk and free its memory
			RC writeStatus = writeBTreeNode(tree, newLeafNode->selfPage, newLeafNode);
			freeBTreeNode(newLeafNode);

			return (status != RC_OK) ? status : writeStatus;
		}
	}
}

// Helper method to insert a key and right child pointer into a non-leaf node
RC insertIntoNode(BTreeHandle* tree, BTreeNode* node, char* serializedKey, int rightChildPageNum) {
	if (tree == NULL || node == NULL || serializedKey == NULL) {
		return RC_FILE_HANDLE_NOT_INIT;
	}

	// This method assumes that the node is not full and has space for the new key
	int i;

	// Find the position to insert the new key
	for (i = node->numKeys; i > 0 && strcmp(node->keyData[i - 1], serializedKey) > 0; i--) {
		node->keyData[i] = node->keyData[i - 1];
		node->children[i + 1] = node->children[i]; // Shift children pointers as well
	}

	// Insert the new key and update the corresponding right child pointer
	node->keyData[i] = strdup(serializedKey); // Duplicate the serialized key
	node->children[i + 1] = rightChildPageNum;
	node->numKeys++;

	// Write the updated node back to the disk
	RC writeStatus = writeBTreeNode(tree, node->selfPage, node);
	return writeStatus;
}

// Helper function to split an overflowing internal node and create a new parent if necessary
RC splitNode(BTreeHandle* tree, BTreeNode* nodeToSplit, BTreeNode** newNode, char* newSerializedKey, int newChildPage) {
	BTreeData* treeData = (BTreeData*)tree->mgmtData;

	// Allocate a new node which will store the keys after the split
	RC rc = createBTreeNode(tree, newNode, false);
	if (rc != RC_OK) return rc;

	// Increment lastPage to get a new page number for the new node
	treeData->lastPage++;
	(*newNode)->selfPage = treeData->lastPage;

	// Find the median key's index
	int medianIndex = nodeToSplit->numKeys / 2;
	bool newKeyRightSide = strcmp(newSerializedKey, nodeToSplit->keyData[medianIndex]) > 0;

	// Determine the split point and which side the new key should be on
	if (newKeyRightSide) {
		medianIndex++;
	}

	// Move half of the keys and children to the new node
	(*newNode)->numKeys = 0;
	int i = medianIndex;
	if (newKeyRightSide) i++;
	for (; i < nodeToSplit->numKeys; i++) {
		(*newNode)->keyData[(*newNode)->numKeys] = strdup(nodeToSplit->keyData[i]); // Copying keys
		(*newNode)->children[(*newNode)->numKeys] = nodeToSplit->children[i];
		(*newNode)->numKeys++;
	}
	(*newNode)->children[(*newNode)->numKeys] = nodeToSplit->children[i];

	// Adjust the number of keys in the original node
	nodeToSplit->numKeys = medianIndex;
	if (newKeyRightSide) nodeToSplit->numKeys--;

	// Insert the new key in the appropriate node
	if (newKeyRightSide) {
		insertIntoNode(tree, *newNode, newSerializedKey, newChildPage);
	}
	else {
		insertIntoNode(tree, nodeToSplit, newSerializedKey, newChildPage);
	}

	// Update children's parentPage
	// ... [Rest of the function remains the same] ...

	return RC_OK;
}

// helper-method for insertKey to split a full leaf node
RC splitLeafNode(BTreeHandle* tree, BTreeNode* leafNode, char* serializedKey, RID rid, char** newKey, int* newPageNum) {
	if (tree == NULL || leafNode == NULL) {
		return RC_FILE_HANDLE_NOT_INIT;
	}

	if (!leafNode->leaf) {
		return RC_IM_KEY_NOT_FOUND;
	}

	BTreeData* treeData = (BTreeData*)tree->mgmtData;
	BTreeNode* newLeafNode;
	RC rc = createBTreeNode(tree, &newLeafNode, true);
	if (rc != RC_OK) {
		return rc;
	}

	// Determine the split point
	int splitPoint = leafNode->numKeys / 2;
	bool insertInNewNode = false;

	// Initialize new leaf node properties
	newLeafNode->numKeys = 0;
	newLeafNode->parentPage = leafNode->parentPage;
	newLeafNode->rightPage = leafNode->rightPage;
	newLeafNode->leftPage = leafNode->selfPage;

	// Increment lastPage to get a new page number for the new node
	treeData->lastPage++;
	newLeafNode->selfPage = treeData->lastPage;
	*newPageNum = newLeafNode->selfPage;

	// Split keys and RIDs between the two nodes and insert the new key
	for (int i = splitPoint; i < leafNode->numKeys; i++) {
		if (!insertInNewNode && strcmp(serializedKey, leafNode->keyData[i]) < 0) {
			// Insert the new key into the new leaf node
			newLeafNode->keyData[newLeafNode->numKeys] = strdup(serializedKey);
			newLeafNode->data[newLeafNode->numKeys] = rid;
			newLeafNode->numKeys++;
			insertInNewNode = true;
		}
		// Move the key from the old node to the new node
		newLeafNode->keyData[newLeafNode->numKeys] = leafNode->keyData[i];
		newLeafNode->data[newLeafNode->numKeys] = leafNode->data[i];
		newLeafNode->numKeys++;
	}

	if (!insertInNewNode) {
		// If the new key was not inserted, insert it now
		insertIntoLeaf(tree, newLeafNode, serializedKey, rid);
	}

	// Update the number of keys in the original leaf node
	leafNode->numKeys = splitPoint;

	// Update sibling pointers for the original leaf node
	leafNode->rightPage = newLeafNode->selfPage;

	// Write the new leaf node to disk
	rc = writeBTreeNode(tree, newLeafNode->selfPage, newLeafNode);
	if (rc != RC_OK) {
		freeBTreeNode(newLeafNode);
		return rc;
	}

	// Write the updated original leaf node to disk
	rc = writeBTreeNode(tree, leafNode->selfPage, leafNode);
	if (rc != RC_OK) {
		freeBTreeNode(newLeafNode);
		return rc;
	}

	// Set the new key for the caller if necessary
	if (newKey != NULL) {
		*newKey = strdup(newLeafNode->keyData[0]);
	}

	freeBTreeNode(newLeafNode);
	return RC_OK;
}

// helper-method for insertKey to insert a new key and page number into the parent node
RC insertIntoParent(BTreeHandle* tree, BTreeNode* childNode, char* newSerializedKey, int newPageNum) {
	if (tree == NULL || childNode == NULL) {
		return RC_FILE_HANDLE_NOT_INIT;
	}

	BTreeData* treeData = (BTreeData*)tree->mgmtData;
	int parentPageNum = childNode->parentPage;

	if (parentPageNum == -1) {
		// Creating a new root if there is no parent
		BTreeNode* newRootNode;
		RC rc = createBTreeNode(tree, &newRootNode, false);
		if (rc != RC_OK) {
			return rc;
		}

		newRootNode->keyData[0] = strdup(newSerializedKey); // Copy the serialized key
		newRootNode->children[0] = childNode->selfPage;
		newRootNode->children[1] = newPageNum;
		newRootNode->numKeys = 1;
		newRootNode->parentPage = -1;

		// Update child nodes' parentPage
		BTreeNode* leftChild, * rightChild;
		readBTreeNode(tree, newRootNode->children[0], &leftChild);
		readBTreeNode(tree, newRootNode->children[1], &rightChild);
		leftChild->parentPage = rightChild->parentPage = newRootNode->selfPage;
		writeBTreeNode(tree, newRootNode->children[0], leftChild);
		writeBTreeNode(tree, newRootNode->children[1], rightChild);
		freeBTreeNode(leftChild);
		freeBTreeNode(rightChild);

		treeData->rootPage = newRootNode->selfPage;
		treeData->treeHeight++;

		rc = writeBTreeNode(tree, newRootNode->selfPage, newRootNode);
		if (rc != RC_OK) {
			freeBTreeNode(newRootNode);
			return rc;
		}

		freeBTreeNode(newRootNode);
	}
	else {
		// Handle inserting into an existing parent
		BTreeNode* parentNode;
		RC rc = readBTreeNode(tree, parentPageNum, &parentNode);
		if (rc != RC_OK) {
			return rc;
		}

		BTreeNode* newParentNode = NULL; // For a new node after split, if needed

		// Update the child's parentPage
		BTreeNode* newChildNode;
		rc = readBTreeNode(tree, newPageNum, &newChildNode);
		if (rc != RC_OK) return rc;
		newChildNode->parentPage = parentNode->selfPage;
		writeBTreeNode(tree, newPageNum, newChildNode);
		freeBTreeNode(newChildNode);

		// Check if insertion or split is required
		if (parentNode->numKeys < treeData->n) {
			rc = insertIntoNode(tree, parentNode, newSerializedKey, newPageNum);
		}
		else {
			rc = splitNode(tree, parentNode, &newParentNode, newSerializedKey, newPageNum);
		}

		if (rc != RC_OK) {
			freeBTreeNode(parentNode);
			if (newParentNode != NULL) {
				freeBTreeNode(newParentNode);
			}
			return rc;
		}

		// Write back the parent node
		rc = writeBTreeNode(tree, parentNode->selfPage, parentNode);
		freeBTreeNode(parentNode);
		if (newParentNode != NULL) {
			freeBTreeNode(newParentNode);
		}
	}

	return RC_OK;
}

extern RC deleteKey(BTreeHandle* tree, Value* key) {
	if (tree == NULL || key == NULL) {
		return RC_IM_KEY_NOT_FOUND;
	}

	// Serialize the key for comparison
	char* serializedKey = serializeValue(key);
	if (!serializedKey) {
		return RC_IM_KEY_NOT_FOUND;  // or a more appropriate error
	}

	BTreeNode* leafNode;

	// Find the leaf node that may contain the key
	RC findStatus = findLeafNode(tree, serializedKey, &leafNode);
	if (findStatus != RC_OK) {
		free(serializedKey);
		return findStatus;
	}

	// Attempt to delete the key from the leaf node
	RC deleteStatus = deleteFromLeaf(leafNode, serializedKey);
	if (deleteStatus != RC_OK) {
		free(serializedKey);
		freeBTreeNode(leafNode);
		return deleteStatus;
	}

	// Handle any underflow that might occur after deletion
	RC underflowStatus = handleUnderflow(tree, leafNode);
	if (underflowStatus != RC_OK) {
		free(serializedKey);
		freeBTreeNode(leafNode);
		return underflowStatus;
	}

	free(serializedKey);
	freeBTreeNode(leafNode);
	return RC_OK;
}
// helper-method to delete a key from a leaf node
RC deleteFromLeaf(BTreeNode* leafNode, char* serializedKey) {
	if (leafNode == NULL) {
		return RC_FILE_HANDLE_NOT_INIT;
	}

	// checking if leaf node is leaf
	if (!leafNode->leaf) {
		return RC_IM_KEY_NOT_FOUND;
	}

	int i;
	for (i = 0; i < leafNode->numKeys; i++) {
		if (strcmp(leafNode->keyData[i], serializedKey) == 0) {
			free(leafNode->keyData[i]);  // Free the memory of the key to be deleted
			break;
		}
	}

	if (i == leafNode->numKeys) {
		return RC_IM_KEY_NOT_FOUND; // Key not found
	}

	// Shift keys and RIDs left to remove the deleted key and its associated RID
	for (; i < leafNode->numKeys - 1; i++) {
		leafNode->keyData[i] = leafNode->keyData[i + 1];
		leafNode->data[i] = leafNode->data[i + 1];
	}

	// Nullify the last key slot to prevent dangling pointer
	leafNode->keyData[leafNode->numKeys - 1] = NULL;

	leafNode->numKeys--;

	return RC_OK;
}
// helper-method to redistribute keys between a leaf node and its sibling
RC redistributeLeafNodes(BTreeHandle* tree, BTreeNode* leafNode, BTreeNode* sibling, bool isLeftSibling) {
	if (tree == NULL || leafNode == NULL || sibling == NULL) {
		return RC_FILE_HANDLE_NOT_INIT;
	}

	BTreeData* treeData = (BTreeData*)tree->mgmtData;
	int totalKeys = leafNode->numKeys + sibling->numKeys;
	int redistributeCount = (totalKeys - leafNode->numKeys) / 2;

	if (isLeftSibling) {
		// Shift keys and RIDs in leafNode to make space for new keys
		memmove(leafNode->keyData + redistributeCount, leafNode->keyData, leafNode->numKeys * sizeof(char*));
		memmove(leafNode->data + redistributeCount, leafNode->data, leafNode->numKeys * sizeof(RID));

		// Copy keys from sibling to leafNode
		for (int i = 0; i < redistributeCount; ++i) {
			leafNode->keyData[i] = strdup(sibling->keyData[sibling->numKeys - redistributeCount + i]);
			leafNode->data[i] = sibling->data[sibling->numKeys - redistributeCount + i];
			free(sibling->keyData[sibling->numKeys - redistributeCount + i]); // Free moved key in sibling
			sibling->keyData[sibling->numKeys - redistributeCount + i] = NULL;
		}

		leafNode->numKeys += redistributeCount;
		sibling->numKeys -= redistributeCount;
	}
	else {
		// Copy keys from sibling to leafNode
		for (int i = 0; i < redistributeCount; ++i) {
			leafNode->keyData[leafNode->numKeys + i] = strdup(sibling->keyData[i]);
			leafNode->data[leafNode->numKeys + i] = sibling->data[i];
			free(sibling->keyData[i]); // Free moved key in sibling
			sibling->keyData[i] = NULL;
		}

		// Shift keys and RIDs in sibling to fill the gap
		memmove(sibling->keyData, sibling->keyData + redistributeCount, (sibling->numKeys - redistributeCount) * sizeof(char*));
		memmove(sibling->data, sibling->data + redistributeCount, (sibling->numKeys - redistributeCount) * sizeof(RID));

		leafNode->numKeys += redistributeCount;
		sibling->numKeys -= redistributeCount;
	}

	// Update the parent key value to reflect changes in its children
	BTreeNode* parentNode;
	RC readStatus = readBTreeNode(tree, leafNode->parentPage, &parentNode);
	if (readStatus != RC_OK) {
		return readStatus;
	}

	// Update the parent key based on whether the sibling is left or right
	char* siblingKey = sibling->keyData[0];
	for (int i = 0; i < parentNode->numKeys; i++) {
		if (parentNode->children[i] == sibling->selfPage) {
			if (!isLeftSibling) {
				free(parentNode->keyData[i - 1]); // Free the old key
				parentNode->keyData[i - 1] = strdup(siblingKey);
			}
		}
		else if (parentNode->children[i] == leafNode->selfPage && isLeftSibling) {
			free(parentNode->keyData[i]); // Free the old key
			parentNode->keyData[i] = strdup(leafNode->keyData[0]);
		}
	}

	RC writeStatus = writeBTreeNode(tree, parentNode->selfPage, parentNode);
	if (writeStatus != RC_OK) {
		freeBTreeNode(parentNode);
		return writeStatus;
	}

	freeBTreeNode(parentNode);
	return RC_OK;
}
// helper-method to handle underflow in a leaf node
RC handleUnderflow(BTreeHandle* tree, BTreeNode* leafNode) {
	if (tree == NULL || leafNode == NULL) {
		return RC_FILE_HANDLE_NOT_INIT;
	}

	BTreeData* treeData = (BTreeData*)tree->mgmtData;
	int minKeys = (treeData->n + 1) / 2;

	if (leafNode->numKeys >= minKeys) {
		return RC_OK;
	}

	BTreeNode* leftSibling = NULL;
	BTreeNode* rightSibling = NULL;

	// Check left sibling
	if (leafNode->leftPage != -1) {
		readBTreeNode(tree, leafNode->leftPage, &leftSibling);
		if (leftSibling != NULL && leftSibling->numKeys > minKeys) {
			redistributeLeafNodes(tree, leafNode, leftSibling, true);
			freeBTreeNode(leftSibling);
			return RC_OK;
		}
	}

	// Check right sibling
	if (leafNode->rightPage != -1) {
		readBTreeNode(tree, leafNode->rightPage, &rightSibling);
		if (rightSibling != NULL && rightSibling->numKeys > minKeys) {
			redistributeLeafNodes(tree, leafNode, rightSibling, false);
			freeBTreeNode(rightSibling);
			return RC_OK;
		}
	}

	// Merge leaf nodes
	if (leftSibling != NULL) {
		mergeLeafNodes(tree, leafNode, leftSibling);
		freeBTreeNode(leftSibling);
	}
	else if (rightSibling != NULL) {
		mergeLeafNodes(tree, leafNode, rightSibling);
		freeBTreeNode(rightSibling);
	}

	// Handle the case where there are no siblings to borrow from
	if (!leftSibling && !rightSibling) {
		if (leafNode->parentPage == -1 && leafNode->numKeys == 0) {
			// Single empty node left in the tree
			return RC_OK;
		}
		return RC_IM_NO_MORE_ENTRIES;
	}

	// Handle underflow in the parent node
	BTreeNode* parentNode;
	RC readStatus = readBTreeNode(tree, leafNode->parentPage, &parentNode);
	if (readStatus != RC_OK) {
		return readStatus;
	}

	RC status = handleUnderflow(tree, parentNode);
	freeBTreeNode(parentNode);
	return status;
}
// helper-method to merge two leaf nodes
RC mergeLeafNodes(BTreeHandle* tree, BTreeNode* leafNode, BTreeNode* sibling) {
	if (tree == NULL || leafNode == NULL || sibling == NULL) {
		return RC_FILE_HANDLE_NOT_INIT;
	}

	BTreeData* treeData = (BTreeData*)tree->mgmtData;
	BTreeNode* parentNode = NULL;

	// Read the parent node
	RC readStatus = readBTreeNode(tree, leafNode->parentPage, &parentNode);
	if (readStatus != RC_OK) {
		return readStatus;
	}

	// Merging keys and RIDs from leafNode to sibling
	for (int i = 0; i < leafNode->numKeys; i++) {
		sibling->keyData[sibling->numKeys + i] = strdup(leafNode->keyData[i]); // Duplicate the key
		sibling->data[sibling->numKeys + i] = leafNode->data[i];
		free(leafNode->keyData[i]); // Free the original key in leafNode
	}
	sibling->numKeys += leafNode->numKeys;

	// Adjust sibling links
	sibling->rightPage = leafNode->rightPage;

	// Remove the reference of leafNode in parentNode
	int index;
	for (index = 0; index < parentNode->numKeys && parentNode->children[index] != leafNode->selfPage; index++);
	for (int i = index; i < parentNode->numKeys - 1; i++) {
		parentNode->keyData[i] = parentNode->keyData[i + 1];
		parentNode->children[i] = parentNode->children[i + 1];
	}
	parentNode->numKeys--;

	// Write the updated sibling node back to the disk
	RC writeStatus = writeBTreeNode(tree, sibling->selfPage, sibling);
	if (writeStatus != RC_OK) {
		freeBTreeNode(parentNode);
		return writeStatus;
	}

	// Write the updated parent node back to the disk
	writeStatus = writeBTreeNode(tree, parentNode->selfPage, parentNode);
	if (writeStatus != RC_OK) {
		freeBTreeNode(parentNode);
		return writeStatus;
	}

	// Free the leafNode and handle underflow in the parent node
	RC deleteStatus = freeBTreeNode(leafNode);
	if (deleteStatus != RC_OK) {
		freeBTreeNode(parentNode);
		return deleteStatus;
	}

	if (parentNode->numKeys < (treeData->n + 1) / 2) {
		RC underflowStatus = handleUnderflow(tree, parentNode);
		if (underflowStatus != RC_OK) {
			freeBTreeNode(parentNode);
			return underflowStatus;
		}
	}

	freeBTreeNode(parentNode);
	return RC_OK;
}
extern RC openTreeScan(BTreeHandle* tree, BT_ScanHandle** handle) {
	if (tree == NULL || handle == NULL) {
		return RC_FILE_HANDLE_NOT_INIT;
	}

	
	*handle = (BT_ScanHandle*)malloc(sizeof(BT_ScanHandle));
	if (*handle == NULL) {
		return RC_FILE_HANDLE_NOT_INIT; 
	}

	
	TreeScanData* scanData = (TreeScanData*)malloc(sizeof(TreeScanData));
	if (scanData == NULL) {
		free(*handle);
		return RC_FILE_HANDLE_NOT_INIT; 
	}

	
	BTreeData* btreeData = (BTreeData*)tree->mgmtData;
	scanData->pageNum = btreeData->rootPage; 
	scanData->index = -1; 
	scanData->BTN = NULL; 

	(*handle)->tree = tree;
	(*handle)->mgmtData = scanData;

	return RC_OK;
}

extern RC nextEntry(BT_ScanHandle* handle, RID* result) {
	if (handle == NULL || result == NULL) {
		return RC_FILE_HANDLE_NOT_INIT;
	}

	TreeScanData* scanData = (TreeScanData*)handle->mgmtData;
	BTreeHandle* tree = handle->tree;
	BTreeData* btreeData = (BTreeData*)tree->mgmtData;

	
	if (scanData->pageNum == -1 || (scanData->BTN != NULL && scanData->index >= scanData->BTN->numKeys - 1)) {
		
		if (scanData->BTN != NULL) {
			freeBTreeNode(scanData->BTN);
			scanData->BTN = NULL;
		}

		
		if (scanData->pageNum == -1) {
			scanData->pageNum = btreeData->rootPage;
			while (TRUE) {
				readBTreeNode(tree, scanData->pageNum, &scanData->BTN);
				if (scanData->BTN->leaf) {
					break;
				}
				else {
					scanData->pageNum = scanData->BTN->children[0];
					freeBTreeNode(scanData->BTN);
				}
			}
		}
		else if (scanData->BTN->rightPage != -1) {
			
			scanData->pageNum = scanData->BTN->rightPage;
			readBTreeNode(tree, scanData->pageNum, &scanData->BTN);
		}
		else {
			
			return RC_IM_NO_MORE_ENTRIES;
		}

		scanData->index = -1; 
	}

	
	scanData->index++;
	*result = scanData->BTN->data[scanData->index];
	return RC_OK;
}

extern RC closeTreeScan(BT_ScanHandle* handle) {
	if (handle == NULL) {
		return RC_FILE_HANDLE_NOT_INIT;
	}

	TreeScanData* scanData = (TreeScanData*)handle->mgmtData;
	if (scanData != NULL) {
		
		if (scanData->BTN != NULL) {
			freeBTreeNode(scanData->BTN);
		}

		
		free(scanData);
	}

	
	free(handle);

	return RC_OK;
}

// debug and test functions
extern char* printTree(BTreeHandle* tree)
{
	if (tree == 0)
	{
		return RC_FILE_NOT_FOUND;
	}
	BTreeData* BTD = (BTreeData*)tree->mgmtData;
	if (BTD->entries == 0)
	{
		return RC_OK;
	}
	BTD->walkTreeIndex = malloc(sizeof(int) * (BTD->lastPage + 1));
	BTD->walkTreeIndex[0] = -1;
	int* index = malloc(sizeof(int));
	index[0] = 0;
	BTD->walkTreeIndex[BTD->rootPage] = index[0];
	index[0]++;
	buildTreeIndex(tree, BTD->rootPage, 0, index);
	printTreeNode(tree, BTD->rootPage);
	printf("Entries: %i\n", BTD->entries);
	free(index);
	free(BTD->walkTreeIndex);
	BTD->walkTreeIndex = 0;
	return RC_OK;
}