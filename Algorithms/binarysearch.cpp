//TEMPLATE FOR BINARY SEARCH
//1.) For int return-type
int binarySearch(vector<int> arr, int target) {
	int array_size = arr.size();
	int low = 0;
	int high = array_size - 1;
	while (low <= high) {
		int mid = (low + (high - low) / 2);
		if (arr[mid] == target) {
			return mid;
		}
		else if (arr[mid] < target) {
			low = mid + 1;
		}
		else {
			high = mid - 1;
		}
	}
	return -1;
}
//2.) For boolean return-type
bool binarySearch() {
	int array_size = arr.size();
	int low = 0;
	int high = array_size - 1;
	while (low <= high) {
		int mid = (low + (high - low) / 2);
		if (arr[mid] == target) {
			return mid;
			//return true;
		}
		else if (arr[mid] < target) {
			low = mid + 1;
		}
		else {
			high = mid - 1;
		}
	}
	return false;
}
//3.) Smaller Code
bool search(int x[], int n, int k) {
	int p = 0;
	for (int a = n; a >= 1; a /= 2) {
		while (p + a < n && x[p + a] <= k) p += a;
	}
	return x[p] == k;
}
// ----------------------------------------------------------------------------------------
/*
COMPLEXITY:- O(log(N))
*/

//Binary Search STL(Standard Template Library)

// function ==> binary_search(startaddress, endaddress, valuetofind);

/*
 int a[] = { 1, 5, 8, 9, 6, 7, 3, 4, 2, 0 };
 if (binary_search(a, a + 10, 10))
    cout << "\nElement found in the array";
 else
    cout << "\nElement not found in the array";

 OUTPUT:- Element not found in the array
*/