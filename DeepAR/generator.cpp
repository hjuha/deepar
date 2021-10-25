#include <bits/stdc++.h>
#include <boost/math/distributions/inverse_gamma.hpp>

#define EPSILON 0.000000001

using namespace std;

template<typename T>
using Matrix = vector<vector<T>>;

mt19937 gen;

template <typename T>
vector<int> hopcroft_karp(Matrix<T> matrix) {
	int n = matrix.size();
	vector<int> bipartite_edges[n];
	int matching[n];
	int inv[n];
	for (int i = 0; i < n; i++) {
		matching[i] = inv[i] = -1;
		bipartite_edges[i].clear();
		for (int j = 0; j < n; j++) {
			if (matrix[i][j] >= EPSILON) {
				bipartite_edges[i].push_back(j + n);
			}
		}
	}

	vector<int> dag[2 * n];
	bool modified = false;
	do {
		modified = false;
		queue<int> q;
		int seen[2 * n];
		for (int i = 0; i < n; i++) {
			seen[i] = seen[i + n] = 0;
			dag[i].clear();
			dag[i + n].clear();
			if (matching[i] == -1) {
				q.push(i);
			}
		}
		// construct graph
		bool stop = false;
		while (!q.empty()) {
			int i = q.front();
			q.pop();
			if (seen[i]) continue;
			seen[i] = 1;
			if (i >= n) {
				if (inv[i - n] == -1) {
					stop = true;
				} else if (!seen[inv[i - n]] && !stop) {
					q.push(inv[i - n]);
					dag[i].push_back(inv[i - n]);
				}
			} else {
				for (int j : bipartite_edges[i]) {
					if (!seen[j] && !stop && matching[i] != j) {
						q.push(j);
						dag[i].push_back(j);
					}
				}
			}
		}
		bool handled[2 * n];
		for (int i = 0; i < 2 * n; i++) {
			seen[i] = -1;
			handled[i] = false;
		}

		stack<int> s;
		for (int i = 0; i < n; i++) {
			if (matching[i] == -1) {
				s.push(i);
				seen[i] = -2;
			}
			while (!s.empty()) {
				int i = s.top();
				s.pop();
				if (handled[i]) continue;
				handled[i] = true;
				if (i >= n && inv[i - n] == -1) {
					while (!s.empty()) {
						s.pop();
					}
					modified = true;
					int j = i;
					while (j != -2) {
						matching[seen[j]] = j;
						inv[j - n] = seen[j];
						j = seen[seen[j]];
					}
					break;
				} else {
					for (int j : dag[i]) {
						if (handled[j]) continue;
						seen[j] = i;
						s.push(j);
					}
				}
			}
		}
	} while (modified);

	vector<int> perm(n);
	for (int i = 0; i < n; i++) {
		perm[i] = matching[i] - n;
	}
	return perm;
}

template<typename T>
vector<int> scc(Matrix<T> m) {
	int n = m.size();
	vector<int> v[n];
	vector<int> r[n];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (m[i][j] > EPSILON) {
				v[i].push_back(j);
				r[j].push_back(i);
			}
		}
	}
	vector<int> state(n);
	vector<int> order;

	for (int i = 0; i < n; i++) {
		if (state[i]) continue;
		stack<int> s;
		s.push(i);
		state[i] = 1;
		while (!s.empty()) {
			int i = s.top();
			if (!v[i].empty()) {
				int j = v[i].back();
				v[i].pop_back();
				if (!state[j]) {
					s.push(j);
					state[j] = 1;
				}
			} else {
				s.pop();
				order.push_back(i);
			}
		}
	}
	vector<int> components(n);
	int component_id = 1;
	reverse(order.begin(), order.end());
	for (int i : order) {
		if (components[i]) continue;
		stack<int> s;
		s.push(i);
		while (!s.empty()) {
			int i = s.top();
			s.pop();
			if (components[i]) {
				continue;
			}
			components[i] = component_id;
			for (int j : r[i]) {
				s.push(j);
			}
		}
		component_id += 1;
	}

	return components;
}

template<typename T>
Matrix<T> tassa(Matrix<T> matrix) {
	int n = matrix.size();
	vector<int> matching = hopcroft_karp(matrix);
	int p[n];
	int pr[n];
	for (int i = 0; i < n; i++) {
		if (matching[i] < 0) return {};
		p[i] = matching[i];
		pr[matching[i]] = i;
		matching[i] = i;
	}
	{
		Matrix<T> permuted;
		for (int i = 0; i < n; i++) {
			permuted.push_back(matrix[pr[i]]);
		}
		matrix = permuted;
	}
	vector<long double> cs(n);
	for (int i = 0; i < n; i++) {
		cs[i] = matrix[i][i];
		matrix[i][i] = 0;
	}
	vector<int> components = scc(matrix);
	Matrix<T> support = matrix;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			support[i][j] = 0;
		}
	}
	for (int i = 0; i < n; i++) {
		support[i][i] = cs[i];
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (matrix[i][j] > EPSILON && components[i] == components[j]) {
				support[i][j] = matrix[i][j];
			}
		}
	}
	{
		Matrix<T> permuted;
		for (int i = 0; i < n; i++) {
			permuted.push_back(support[p[i]]);
		}
		support = permuted;
	}
	return support;
}

int main() {
	int n;
	long double p;
	cin>>n>>p;
	cout<<n<<endl;
	srand((int)(n + 10000 / p));
	while (true) {
		Matrix<long double> matrix;
		for (int i = 0; i < n; i++) {
			matrix.push_back(vector<long double>(n));
			for (int j = 0; j < n; j++) {
				if ((rand() % 1000000) / 1000000.0 < p) matrix[i][j] = 1;
			}
		}
		if (tassa(matrix).empty()) continue;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cout<<matrix[i][j];
				if (j != n - 1) cout<<" ";
			}
			cout<<endl;
		}
		break;
	}
}