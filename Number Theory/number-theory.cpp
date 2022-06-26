/* 1.) GCD OR HCF (GREATEST COMMON DIVISOR) */
#include <bits/stdc++.h>
using namespace std;

int gcd(int a, int b) {
	if (a == 0) return b;
	return (b % a, a);
}
//Complexity:- O(log(min(a,b)))

int gcd_stl = __gcd(a, b); // Complexity:- O(log(2n))

int gcd_Extended(int a, int b, int *x, int *y) {
	if (a == 0) {
		*x = 0;
		*y = 1;
		return b;
	}
	int x1, y1;
	int gcd = gcdExtended(b % a, a, &x1, &y1);
	*x = y1 - (b / a) * x1;
	*y = x1;
	return gcd;
}
//Complexity:- O(log(max(a,b)))

// //----------------------------------------------------------------------------------------------------------------------------------
/*LCM or (Least Common Factor)*/

int a, b; cin >> a >> b;
int lcm = (a*b) / gcd(a, b);
//Complexity:- O(log(max(a,b)))
//Auxillary Space:- O(log(max(a,b)))

//----------------------------------------------------------------------------------------------------------------------------------

/*FACTORIAL OF A NUMBER*/

int factorial(int n) {
	return (n == 1 || n == 0) ? 1 : n * factorial(n - 1);
}
//Complexity:- O(n)

//LARGE NUMBER FACTORIAL
int multiply(int x, int res[], int res_size);

void factorial(int n) {
	int res[MAX];
	res[0] = 1;
	int res_size = 1;
	for (int x = 2; x <= n; x++) {
		res_size = multiply(x, res, res_size);
	}
	cout << "Factorial of given number is \n";
	for (int i = res_size - 1; i >= 0; i--)
		cout << res[i];
}

int multiply(int x, int res[], int res_size) {
	int carry = 0;
	for (int i = 0; i < res_size; i++) {
		int prod = res[i] * x + carry;
		res[i] = prod % 10;
		carry  = prod / 10;
	}
	while (carry) {
		res[res_size] = carry % 10;
		carry = carry / 10;
		res_size++;
	}
	return res_size;
}
//Complexity:- O(nlog(n!))

//----------------------------------------------------------------------------------------------------------------------------------

// SIEVE OF ErRATOSTHENES AND SEGMENTED SIEVE

//First Approach (Normally)
void sieve() {
	bool prime[n + 1];
	memset(prime, true, sizeof(prime));
	for (int i = 2; i * i <= n; ++i) {
		if (prime[p] == true) {
			for (int j = i * i; j <= n; j += i) {
				prime[i] = false;
			}
		}
	}
	for (int i = 2; i <= n; ++i) {
		if (prime[p]) {
			cout << p << " ";
		}
	}
}

//Second Approach (BitSet method)
#include <bitset> //this header should be included especially if you are not using  #include <bits/stdc++.h>

bitset<500001> prime;
void sieve(int n) {
	prime[0] = 1;
	for (int i = 2; i <= n; i += 2) {
		if (prime[i / 2] == 0) {
			for (int i = 3 * i; j <= n; j += 2 * i) {
				prime[j / 2] = 1;
			}
		}
	}
}

//Time Complexity:- O(n*log(log(n)))

//SEGMENTED SIEVE

void sieve_normal(int limit, vecotr<int> &prime) {
	vector<bool> mark(limit + 1, true);
	for (int i = 2; i * i < limit; i++) {
		if (mark[i] == true) {
			for (int j = i * i; j < limit; j += p) {
				mark[i] = false;
			}
		}
	}
	for (int i = 2; i < n; i++) {
		if (mark[i] == true) {
			prime.push_back(i);
			cout << p << " ";
		}
	}
}

void segmented_sieve() {
	int limit = floor(sqrt(n)) + 1;
	vector<int> prime;
	prime.reserve(limit);
	simpleSieve(limit, prime);
	int low = limit;
	int high = 2 * limit;
	while (low < n) {
		if (high >= n) {
			high = n;
		}
		bool mark[limit + 1];
		memset(mark, true, sizeof(mark)) {
			for (int i = 0; i < prime.size(); i++) {
				int loLim = floor(low / prime[i]) * prime[i];
				if (loLim < low) {
					loLim += prime[i];
				}
			}{}
			for (int j = loLim; j < high; j += prime[i]) {
				mark[j - low] = false;
			}
			for (int i = low; i < high; i++)
				if (mark[i - low] == true)
					cout << i << " ";
			low = low + limit;
			high = high + limit;
		}
	}

//----------------------------------------------------------------------------------------------------------------------------------

//Effective way for finding area of an individual rectangle is W*L

	long long area(int bl_x, int bl_y, int tr_x, tr_y) {
		long long length = (tr_y - bl_y);
		long long width  = (tr_x - bl_x);
		return (length * width);
	}


//CHECKING IF TWO RECTANGLES INTERSECT EACH OTHER

// Given two rectangles a and b, there are only two cases where they do not intersect:
//tr_a_y <= bl_b_y OR bl_a_y >= tr_b_y
//br_a_x >= tr_b_x OR tr_a_x <= bl_b_x

	vector<int> intersect(vector<int> s1, vector<int> s2) {
		int bl_a_x = s1[0], bl_a_y = s1[1], tr_a_x = s1[2], tr_a_y = s1[3];
		int bl_b_x = s2[0], bl_b_y = s2[1], tr_b_x = s2[2], tr_b_y = s2[3];

		// no overlap
		if (bl_a_x >= tr_b_x || tr_a_x <= bl_b_x
		        || bl_a_y >= tr_b_y || tr_a_y <= bl_b_y) {
			return false;
		} else {
			return true;
		}
	}

//FINDING AREA OF INTERSECTION

//width    = min(tr_a_x, tr_b_x) - max(bl_a_x,bl_b_x)
//length = min(tr_a_y, tr_b_y) - max(bl_a_y,bl_b_y);

	int inter_area(vector<int> s1, vector<int> s2) {
		int bl_a_x = s1[0], bl_a_y = s1[1], tr_a_x = s1[2], tr_a_y = s1[3];
		int bl_b_x = s2[0], bl_b_y = s2[1], tr_b_x = s2[2], tr_b_y = s2[3];

		return (
		           (min(tr_a_x, tr_b_x) - max(bl_a_x, bl_b_x))
		           * (min(tr_a_y, tr_b_y) - max(bl_a_y, bl_b_y))
		       );
	}


//----------------------------------------------------------------------------------------------------------------------------------

//BINOMIAL CO-EFFICIENT

//C(n, k) ==> coefficient of x^k in the expansion of (1 + x)^n.

	/*
	1.) Optimal Sub-Structure
	C(n, k) = C(n-1, k-1) + C(n-1, k)
	C(n, 0) = C(n, n) = 1
	*/

	int binomialCoeff(int n, int k)
	{
		// Base Cases
		if (k > n)
			return 0;
		if (k == 0 || k == n)
			return 1;

		return (binomialCoeff(n - 1, k - 1) + binomialCoeff(n - 1, k));
	}

//2.) 2nd Approach
	int binomialCoeff(int n, int k)
	{
		int C[n + 1][k + 1];
		int i, j;

		// Calculate value of Binomial Coefficient
		// in bottom up manner
		for (i = 0; i <= n; i++) {
			for (j = 0; j <= min(i, k); j++) {
				// Base Cases
				if (j == 0 || j == i)
					C[i][j] = 1;

				// Calculate value using previously
				// stored values
				else
					C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
			}
		}

		return C[n][k];
	}
// Time Complexity and Space Complexity:- O(n*k)

//3.) 3rd Approach
	int binomialCoeff(int n, int k)
	{
		int C[k + 1];
		memset(C, 0, sizeof(C));

		C[0] = 1; // nC0 is 1

		for (int i = 1; i <= n; i++)
		{

			// Compute next row of pascal triangle using
			// the previous row
			for (int j = min(i, k); j > 0; j--)
				C[j] = C[j] + C[j - 1];
		}
		return C[k];
	}
// Time Complexity:- O(n*k) and Space Complexity:- O(k)





//By GCD

	int gcd(int a, int b)
	{
		if (b == 0)
			return a;
		return gcd(b, a % b);
	}

	int nCr(int n, int r)
	{
		// base case
		if (r > n)
			return 0;

		// C(n,r) = C(n,n-r)
		if (r > n - r)
			r = n - r;

		int mod = 1000000007;

		// array of elements from n-r+1 to n
		int arr[r];

		for (int i = n - r + 1; i <= n; i++) {
			arr[i + r - n - 1] = i;
		}

		long int ans = 1;
		// for numbers from 1 to r find arr[j],
		// such that gcd(i,arr[j])>1
		for (int k = 1; k < r + 1; k++) {
			int j = 0, i = k;
			while (j < r) {
				int x = gcd(i, arr[j]);
				if (x > 1) {
					// if gcd>1, divide both by gcd
					arr[j] /= x;
					i /= x;
				}

				// if i becomes 1, no need to search arr
				if (i == 1)
					break;
				j += 1;
			}
		}

		// single pass to multiply the numerator
		for (int i : arr)
			ans = (ans * i) % mod;
		return (int)ans;
	}

//Time Complexity:-  O(( min(r, n-r)^2 ) * log(n)) and Space Complexity:- O(min(r, n-r))

//----------------------------------------------------------------------------------------------------------------------------------

//Finding remainder without modulo operator

	int getRemainder(int num, int divisor) {
		return (num - divisor * (num / divisor));
	}


//----------------------------------------------------------------------------------------------------------------------------------

//Modulo Multiplicative Inverse

//	a x ≅ 1 (mod m)

	int modInverse(int a, int m) {
		for (int x = 1; x < m; x++) {
			if (((a % m) * (x % m)) % m == 1) {
				return x;
			}
		}
	}
// Time Complexity:- O(m)


// Method 2 :- Fermat’s Little theorem, O(Log m) [Works when ‘m’ is prime]

	void modInverse(int a, int m) {
		int g = __gcd(a, m);
		if (g != 1)
			cout << "Inverse doesn't exist";
		else
		{
			cout << "Modular multiplicative inverse is " << power(a, m - 2, m);
		}
	}
	int power(int x, unsigned int y, unsigned int m) {
		if (y == 0)
			return 1;
		int p = power(x, y / 2, m) % m;
		p = (p * p) % m;

		return (y % 2 == 0) ? p : (x * p) % m;
	}

//Time Complexity: O(Log m)

//----------------------------------------------------------------------------------------------------------------------------------

// Compute n! under modulo p

	int modFact(int n, int p) {
		if (n >= p)
			return 0;

		int result = 1;
		for (int i = 1; i <= n; i++)
			result = (result * i) % p;

		return result;
	}
//Time Complexity: O(n)

//----------------------------------------------------------------------------------------------------------------------------------

int main(){
	return 0;
}