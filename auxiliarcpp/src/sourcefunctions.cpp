//#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <queue>
//#include <Eigen>;
#include <random>
#include <RcppArmadillo.h>
#include <cmath>    // For pow()
#include <tuple>
#include <unordered_set>
#include <limits>
#include <utility> // For std::pair
//#include <stan/math.hpp>
using namespace Rcpp;
using namespace std;
using namespace arma;

//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
 // [[Rcpp::export]]
 Rcpp::NumericMatrix primSpanningTree(int n) {
   // Initialize adjacency matrix for the spanning tree (n x n)
   Rcpp::NumericMatrix adj_matrix(n, n);
   adj_matrix.fill(0);  // Ensure all values are initialized to zero

   // Random number generator
   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_real_distribution<> dis(0.0, 1.0);

   // Priority queue to store edges (weight, (u, v))
   using Edge = std::pair<double, std::pair<int, int>>;
   std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;

   // Start with node 0 in the tree
   std::vector<bool> in_tree(n, false);
   in_tree[0] = true;

   // Add all edges from node 0 to the priority queue
   for (int v = 1; v < n; ++v) {
     double weight = dis(gen);  // Random weight for the edge
     pq.push({weight, {0, v}});
   }

   // While there are fewer than n-1 edges in the tree
   int edge_count = 0;
   while (edge_count < n - 1 && !pq.empty()) {
     // Get the edge with the smallest weight
     Edge edge = pq.top();
     pq.pop();

     int u = edge.second.first;
     int v = edge.second.second;

     // If both vertices are already in the tree, skip (it would form a cycle)
     if (in_tree[u] && in_tree[v]) {
       continue;
     }

     // Add edge to the tree
     adj_matrix(u, v) = 1;  // Add the edge from u to v
     adj_matrix(v, u) = 1;  // Add the edge from v to u (undirected graph)
     edge_count++;

     // Add the new vertex to the tree and mark it
     int new_vertex = in_tree[u] ? v : u;
     in_tree[new_vertex] = true;

     // Add edges from this new vertex to the priority queue
     for (int w = 0; w < n; ++w) {
       // Check if there is an edge and w is not already in the tree
       if (!in_tree[w]) {
         double new_weight = dis(gen);  // Random weight for the edge
         pq.push({new_weight, {new_vertex, w}});
       }
     }
   }

   return adj_matrix;
 }


/* Build W from adjacency (dense) */
static arma::mat build_W_dense(const arma::mat& A, bool row_standardize) {
  arma::mat B = A;
  B.diag().zeros();
  B = 0.5 * (B + B.t());
  // binarize
  B.for_each([](double& v){ v = (v > 0.0) ? 1.0 : 0.0; });
  if (row_standardize) {
    arma::vec rs = arma::sum(B, 1);
    for (arma::uword i = 0; i < B.n_rows; ++i) {
      if (rs(i) > 0.0) B.row(i) /= rs(i);
    }
  }
  return B;
}

//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
bool is_spanning_tree_cpp(IntegerMatrix adj_matrix) {
  int n = adj_matrix.nrow();  // Number of nodes (rows in the adjacency matrix)
  int edge_count = 0;

  // Use a vector to track visited nodes and simulate a stack
  std::vector<bool> visited(n, false);
  std::vector<int> stack;

  // Start DFS from the first node
  stack.push_back(0);
  visited[0] = true;
  int visited_count = 1;

  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();

    // Loop through neighbors
    for (int neighbor = 0; neighbor < n; ++neighbor) {
      if (adj_matrix(node, neighbor) == 1) {
        ++edge_count; // Count every edge

        if (!visited[neighbor]) {
          visited[neighbor] = true;
          stack.push_back(neighbor);
          ++visited_count;
        }
      }
    }
  }

  // Since each edge is counted twice in the adjacency matrix
  edge_count /= 2;

  // Check if the graph is connected and acyclic (n-1 edges)
  return (visited_count == n && edge_count == (n - 1));
}



//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::depends(RcppArmadillo)]]
inline bool ar2_stationary(const double phi1, const double phi2) {
  // Stationarity triangle:
  // phi2 in (-1,1), phi1+phi2 < 1, phi2 - phi1 < 1
  return (phi2 > -1.0 && phi2 < 1.0 &&
          (phi1 + phi2) < 1.0 &&
          (phi2 - phi1) < 1.0);
}

// [[Rcpp::export]]
arma::mat corr_ar2_fast(const int n, const double phi1, const double phi2,
                        const bool use_clamp = true) {
  if (n < 1) Rcpp::stop("n must be >= 1");
  if (!ar2_stationary(phi1, phi2)) {
    Rcpp::stop("Non-stationary AR(2): need phi2 in (-1,1), phi1+phi2 < 1, phi2-phi1 < 1.");
  }

  // Build autocorrelation vector rho[0..n-1]
  arma::vec rho(n, arma::fill::zeros);
  rho[0] = 1.0;
  if (n > 1) {
    double r1 = phi1 / (1.0 - phi2);     // Yule–Walker
    if (use_clamp) {
      if (r1 >  1.0) r1 =  1.0 - 1e-12;
      if (r1 < -1.0) r1 = -1.0 + 1e-12;
    }
    rho[1] = r1;
  }
  for (int k = 2; k < n; ++k) {
    double v = phi1 * rho[k - 1] + phi2 * rho[k - 2];
    if (use_clamp) {
      if (v >  1.0) v =  1.0 - 1e-12;
      if (v < -1.0) v = -1.0 + 1e-12;
    }
    rho[k] = v;
  }

  // Allocate result
  arma::mat C(n, n, arma::fill::ones);

  // Fill by columns: column j has:
  //  - lower part (i >= j): C(i,j) = rho[i-j]  -> contiguous block copy
  //  - upper part (i <  j): C(i,j) = rho[j-i]  -> reverse copy of rho[1..j]
  // We then mirror to C(j,i) implicitly by writing both halves in-place.

  // Parallelize over columns, if OpenMP available
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j = 0; j < n; ++j) {
    // Lower triangle including diagonal: i = j..n-1
    // length = n - j, copy rho[0..n-j-1] into C(j..n-1, j)
    double* colPtr = C.colptr(j);  // start of column j (contiguous)
    std::memcpy(colPtr + j, rho.memptr(), (size_t)(n - j) * sizeof(double));

    // Upper triangle: i = 0..j-1, need rho[j-i] (reverse of rho[1..j])
    // Write sequentially to keep cache-friendly pattern on the column
    for (int i = 0; i < j; ++i) {
      colPtr[i] = rho[j - i];
    }
  }

  // Ensure symmetry explicitly (numerical safety, very cheap)
  // Only needed if you want to be ultra-safe; you can omit it.
  // C = arma::symmatu(C);

  return C;
}









//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
bool is_connected(const NumericMatrix& adj_matrix) {
  int n = adj_matrix.nrow();
  vector<bool> visited(n, false);
  queue<int> q;

  // Start BFS from node 0
  q.push(0);
  visited[0] = true;
  int visited_count = 1;

  while (!q.empty()) {
    int node = q.front();
    q.pop();

    for (int i = 0; i < n; ++i) {
      if (adj_matrix(node, i) == 1 && !visited[i]) {
        visited[i] = true;
        q.push(i);
        visited_count++;
      }
    }
  }

  return visited_count == n;  // If all nodes are visited, the graph is connected
}

// Helper function to shuffle and select an edge randomly with probability p
pair<int, int> select_edge(vector<pair<int, int>>& possible_edges, double p) {
  // Shuffle the edges to introduce randomness
  random_shuffle(possible_edges.begin(), possible_edges.end());

  // Select an edge based on probability p
  for (const auto& edge : possible_edges) {
    if (R::runif(0, 1) <= p) {
      return edge;
    }
  }

  // If no edge is selected, return a random edge
  return possible_edges[R::runif(0, possible_edges.size())];
}

//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericMatrix random_spanning_tree_prob_2(NumericMatrix adj_matrix, double p) {
  int n = adj_matrix.nrow();

  // Check if the graph is connected
  if (!is_connected(adj_matrix)) {
    stop("The graph is not connected.");
  }

  // Initialize the spanning tree as an empty adjacency matrix
  NumericMatrix spanning_tree(n, n);

  // Set of nodes to visit (starting with node 1, which is index 0 in C++)
  vector<bool> in_tree(n, false);
  in_tree[0] = true;
  int nodes_in_tree = 1;

  // Continue until all nodes are in the spanning tree
  while (nodes_in_tree < n) {
    // Get all edges connecting the current tree to other nodes
    vector<pair<int, int>> possible_edges;
    for (int u = 0; u < n; ++u) {
      if (in_tree[u]) {
        for (int v = 0; v < n; ++v) {
          if (adj_matrix(u, v) == 1 && !in_tree[v]) {
            possible_edges.push_back(make_pair(u, v));
          }
        }
      }
    }

    // Select an edge with probability p
    pair<int, int> selected_edge = select_edge(possible_edges, p);

    // Add the selected edge to the spanning tree
    int u = selected_edge.first;
    int v = selected_edge.second;
    spanning_tree(u, v) = 1;
    spanning_tree(v, u) = 1;

    // Add the new node to the nodes_in_tree set
    in_tree[v] = true;
    nodes_in_tree++;
  }

  // Ensure we have exactly n-1 edges
  int edge_count = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (spanning_tree(i, j) == 1) {
        edge_count++;
      }
    }
  }

  if (edge_count != n - 1) {
    stop("Spanning tree generation failed to connect all nodes.");
  }

  return spanning_tree;
}
// Declare rowsums_mine function
arma::vec rowsums_mine(const arma::mat& a);
//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::vec rowsums_mine(const arma::mat& a) {
  return arma::sum(a, 1);  // Row-wise sum
}
// Declare rowsums_mine function
arma::mat matrixb2_cpp(arma::mat adjmatinf, double rho);
//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat matrixb2_cpp(arma::mat adjmatinf, double rho) {
  arma::vec nneigh = rowsums_mine(adjmatinf);
  arma::vec b = rho / (1 + nneigh * pow(rho, 2));

  for (size_t i = 0; i < nneigh.n_elem; ++i) {
    // Find indices where adjmatinf is 1
    arma::uvec indices = find(adjmatinf.row(i) == 1);
    // Replace those indices with the corresponding value from b
    for (size_t j = 0; j < indices.n_elem; ++j) {
      adjmatinf(i, indices(j)) = b(i); // Directly assign the value from b
    }
  }

  return adjmatinf;
}


arma::mat ar1_correlation_matrix_2(double phi, int lag);
//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat ar1_correlation_matrix_2(double phi, int lag) {
  // int n = lags.size(); // Get the number of lags
  arma::mat corr_matrix(lag, lag, fill::zeros); // Initialize an n x n matrix filled with zeros

  // Fill the correlation matrix using AR(1) correlation formula
  for (int i = 0; i < lag; ++i) {
    for (int j = 0; j < lag; ++j) {
      corr_matrix(i, j) = pow(phi, abs( i- j));
    }
  }

  return corr_matrix; // Return the computed correlation matrix
}





arma::mat matrixF_cpp(const arma::vec& nneigh, double rho);
//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat matrixF_cpp(const arma::vec& nneigh, double rho) {
  arma::vec tau = (1 + (nneigh - 1) * pow(rho, 2)) / (1 - pow(rho, 2));
  return arma::diagmat(tau); // Create a diagonal matrix
}

arma::mat matrixFinv_cpp(const arma::vec& nneigh, double rho);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat matrixFinv_cpp(const arma::vec& nneigh, double rho) {
  arma::vec tau = (1 + (nneigh - 1) * pow(rho, 2)) / (1 - pow(rho, 2));
  return arma::diagmat(1/tau); // Create a diagonal matrix
}

arma::mat varcovspadagar_cpp(const arma::mat& adjmatinf, double rho);
//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat varcovspadagar_cpp(const arma::mat& adjmatinf, double rho) {
  arma::vec nneigh = rowsums_mine(adjmatinf);
  arma::mat matB = matrixb2_cpp(adjmatinf, rho);
  arma::mat matFinv = matrixFinv_cpp(nneigh, rho);

  arma::mat I = arma::eye<arma::mat>(matFinv.n_rows, matFinv.n_cols); // Identity matrix
  arma::mat I_minus_B = I - matB; // Cache I - matB
  arma::mat inv_I_minus_B = arma::inv(I_minus_B);
  arma::mat Sigma = inv_I_minus_B * matFinv * inv_I_minus_B.t();
 // double standard = 1/ trace(Sigma);
  // arma::vec diag_Sigma = Sigma.diag();
  //arma::vec inv_sqrt_diag = 1.0 / sqrt(diag_Sigma);
  //arma::mat D = inv_sqrt_diag * inv_sqrt_diag.t();
 // return Sigma % D;  //+ tau2 * I; // Matrix multiplication to return the result
//return standard*Sigma;
 return Sigma;
}


double trvarcovspadagar_cpp(const arma::mat& adjmatinf, double rho);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
// [[Rcpp::export]]
double trvarcovspadagar_cpp(const arma::mat& adjmatinf, double rho) {
  arma::vec nneigh = rowsums_mine(adjmatinf);
  arma::mat matB = matrixb2_cpp(adjmatinf, rho);
  arma::mat matFinv = matrixFinv_cpp(nneigh, rho);

  arma::mat I = arma::eye<arma::mat>(matFinv.n_rows, matFinv.n_cols); // Identity matrix
  arma::mat I_minus_B = I - matB; // Cache I - matB
  arma::mat inv_I_minus_B = arma::inv(I_minus_B);
  arma::mat Sigma = inv_I_minus_B * matFinv * inv_I_minus_B.t();
  double standard = 1/ trace(Sigma);
  // arma::vec diag_Sigma = Sigma.diag();
  //arma::vec inv_sqrt_diag = 1.0 / sqrt(diag_Sigma);
  //arma::mat D = inv_sqrt_diag * inv_sqrt_diag.t();
  // return Sigma % D;  //+ tau2 * I; // Matrix multiplication to return the result
  return standard;
}




arma::mat varcovspacar_cpp(const arma::mat& adjmatinf, double rho);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::mat varcovspacar_cpp(const arma::mat& adjmatinf, double rho) {
   arma::vec nneigh = rowsums_mine(adjmatinf);
   arma:: mat Nneigh = arma::diagmat(nneigh);
   arma::mat Q = Nneigh - rho*adjmatinf;
   arma::mat invQ = arma::inv(Q);
   arma::vec diag_invQ = invQ.diag();
   arma::vec inv_sqrt_diag = 1.0 / sqrt(diag_invQ);
   arma::mat D = inv_sqrt_diag * inv_sqrt_diag.t();
   return invQ % D;  //+ tau2 * I; // Matrix multiplication to return the result
   //return diag_invQ;
 }



arma::mat varcovdagar_cpp(const arma::mat& adjmatinf, double rho, double sigma2, double tau2);
//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat varcovdagar_cpp(const arma::mat& adjmatinf, double rho, double sigma2, double tau2) {
  arma::vec nneigh = rowsums_mine(adjmatinf);
  arma::mat matB = matrixb2_cpp(adjmatinf, rho);
  arma::mat matF = matrixF_cpp(nneigh, rho);

  arma::mat I = arma::eye<arma::mat>(matF.n_rows, matF.n_cols); // Identity matrix

  return sigma2 * (I - matB).t() * matF * (I - matB) + tau2 * I; // Covariance matrix computation
}


arma::mat vdagar_cpp(const arma::mat& adjmatinf, double rho, double psi);
//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat vdagar_cpp(const arma::mat& adjmatinf, double rho, double psi) {
  arma::vec nneigh = rowsums_mine(adjmatinf);
  arma::mat matB = matrixb2_cpp(adjmatinf, rho);
  arma::mat matFinv = matrixFinv_cpp(nneigh, rho);

  arma::mat I = arma::eye<arma::mat>(matFinv.n_rows, matFinv.n_cols); // Identity matrix
  arma::mat I_minus_B = I - matB; // Cache I - matB
  arma::mat inv_I_minus_B = arma::inv(I_minus_B);

  return inv_I_minus_B * matFinv * inv_I_minus_B.t() + psi * I; // Matrix multiplication to return the result
}

arma::mat spatimecovar_2_cpp(int lag, const arma::mat& adjmatinf,double rho,double psi,
                             double phi);
//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat spatimecovar_2_cpp(int lag, const arma::mat& adjmatinf,double rho,double psi,
                             double phi){

  arma::mat A = varcovspadagar_cpp(adjmatinf,rho);
  arma::mat B =  ar1_correlation_matrix_2(phi,lag);
  arma::mat timespa = arma::kron(A,B);
  arma::mat I = arma::eye<arma::mat>(timespa.n_rows, timespa.n_cols);
  return timespa + psi*I;
}


//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::mat sar_covariance_rcpp(const arma::mat& A_or_W,
                               const double rho,
                               const bool is_W = false,
                               const bool row_standardize = false,
                               const bool return_precision = false) {
   const arma::uword n = A_or_W.n_rows;
   if (A_or_W.n_cols != n) stop("Input must be square.");

   arma::mat W = is_W ? A_or_W : build_W_dense(A_or_W, row_standardize);

   arma::mat I = arma::eye<arma::mat>(n, n);
   arma::mat A1 = I - rho * W;          // (I - rho W)
   arma::mat Q  = A1.t() * A1;          // precision, SPD

   if (return_precision) return Q;

   // Covariance = inv_sympd(Q)
   arma::mat Sigma = arma::inv_sympd(Q);
   //arma::mat Sigma1 = Sigma/arma::trace(Sigma);
   return Sigma;
 }


arma::mat spatimecovarsar_2_cpp(int lag, const arma::mat& adjmat,double rho,double psi,
                                double phi);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::mat spatimecovarsar_2_cpp(int lag, const arma::mat& adjmat,double rho,double psi,
                              double phi){

   arma::mat A = sar_covariance_rcpp(adjmat,rho,false,true,false);
   arma::mat B =  ar1_correlation_matrix_2(phi,lag);
   arma::mat timespa = arma::kron(A,B);
   arma::mat I = arma::eye<arma::mat>(timespa.n_rows, timespa.n_cols);
  return timespa + psi*I;

 }



Rcpp::List spatimecovarsar_211_cpp(int lag, const arma::mat& adjmat,double rho,double psi,
                                double phi);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
Rcpp::List spatimecovarsar_211_cpp(int lag, const arma::mat& adjmat,double rho,double psi,
                                 double phi){

   arma::mat A = sar_covariance_rcpp(adjmat,rho,false,true,false);
   arma::mat B =  ar1_correlation_matrix_2(phi,lag);
   arma::mat timespa = arma::kron(A,B);
   arma::mat I = arma::eye<arma::mat>(timespa.n_rows, timespa.n_cols);
   return Rcpp::List::create(
     Rcpp::Named("A") = A,
     Rcpp::Named("B") = B,
     Rcpp::Named("timespa") = timespa
   );

 }

arma::mat spatimecovarsar_21_cpp(int lag, const arma::mat& adjmatinf,double rho,double psi,
                                double phi);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::mat spatimecovarsar_21_cpp(int lag, const arma::mat& adjmatinf,double rho,double psi,
                                 double phi){

   arma::mat A = sar_covariance_rcpp(adjmatinf,rho,false,true,false);
   arma::mat B =  ar1_correlation_matrix_2(phi,lag);
   arma::mat timespa = arma::kron(A,B);
   arma::mat I = arma::eye<arma::mat>(timespa.n_rows, timespa.n_cols);
   return timespa + psi*I;
 }






arma::mat spatimecovar_ar2_cpp(int lag, const arma::mat& adjmatinf,double rho,double psi,
                             double phi1, double phi2);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::mat spatimecovar_ar2_cpp(int lag, const arma::mat& adjmatinf,double rho,double psi,
                              double phi1, double phi2){

   arma::mat A = varcovspadagar_cpp(adjmatinf,rho);
   arma::mat B =  corr_ar2_fast(lag,phi1,phi2,true);
   arma::mat timespa = arma::kron(A,B);
   arma::mat I = arma::eye<arma::mat>(timespa.n_rows, timespa.n_cols);
   return timespa + psi*I;
 }


arma::mat spatimecovarsar_ar2_cpp(int lag, const arma::mat& adjmat,double rho,double psi,
                               double phi1, double phi2);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::mat spatimecovarsar_ar2_cpp(int lag, const arma::mat& adjmat,double rho,double psi,
                                double phi1, double phi2){

   arma::mat A = sar_covariance_rcpp(adjmat,rho,false,true,false);
   arma::mat B =  corr_ar2_fast(lag,phi1,phi2,true);
   arma::mat timespa = arma::kron(A,B);
   arma::mat I = arma::eye<arma::mat>(timespa.n_rows, timespa.n_cols);
   return timespa + psi*I;
 }


arma::mat spatimecovarsar_2_ar2_cpp(int lag, const arma::mat& adjmatinf,double rho,double psi,
                                  double phi1, double phi2);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
arma::mat spatimecovarsar_2_ar2_cpp(int lag, const arma::mat& adjmatinf,double rho,double psi,
                                   double phi1, double phi2){

   arma::mat A = sar_covariance_rcpp(adjmatinf,rho,false,true,false);
   arma::mat B =  corr_ar2_fast(lag,phi1,phi2,true);
   arma::mat timespa = arma::kron(A,B);
   arma::mat I = arma::eye<arma::mat>(timespa.n_rows, timespa.n_cols);
   return timespa + psi*I;
 }






arma::mat spatimecovarcar_2_cpp(int lag, const arma::mat& adjmatinf,double rho,double psi,
                             double phi);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::mat spatimecovarcar_2_cpp(int lag, const arma::mat& adjmatinf,double rho,double psi,
                              double phi){

   arma::mat A = varcovspacar_cpp(adjmatinf,rho);
   arma::mat B =  ar1_correlation_matrix_2(phi,lag);
   arma::mat timespa = arma::kron(A,B);
   arma::mat I = arma::eye<arma::mat>(timespa.n_rows, timespa.n_cols);
   return timespa + psi*I;
 }




arma::mat spatimecovar_2_zero_cpp(int lag, const arma::mat& adjmatinf,double rho,double psi,
                             double phi);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::mat spatimecovar_2_zero_cpp(int lag, const arma::mat& adjmatinf,double rho,
                              double phi){

   arma::mat A = varcovspadagar_cpp(adjmatinf,rho);
   arma::mat B =  ar1_correlation_matrix_2(phi,lag);
   arma::mat timespa = arma::kron(A,B);
   return timespa;
 }





arma::mat spatimecovar_cpp(int lag, const arma::mat& adjmatinf,double sigma2,double tau2,
                           double phi, double rho);

//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat spatimecovar_cpp(int lag, const arma::mat& adjmatinf,double sigma2,double tau2,
                           double phi, double rho){

  arma::mat A =  ar1_correlation_matrix_2(phi,lag);
  arma::mat B = varcovspadagar_cpp(adjmatinf,rho);
  arma::mat timespa = arma::kron(A,B);
  arma::mat I = arma::eye<arma::mat>(timespa.n_rows, timespa.n_cols);
  return sigma2*timespa + tau2*I;
}

//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double dbeta_rep_cpp(double x, double mu, double sigma2) {
  // Ensure mu is between 0 and 1 to avoid invalid beta distribution parameters
  if (mu <= 0 || mu >= 1) {
    Rcpp::stop("mu must be in (0, 1)");
  }
  // Ensure sigma2 is positive and less than mu(1-mu)
  if (sigma2 <= 0 || sigma2 >= mu * (1 - mu)) {
    Rcpp::stop("Invalid sigma2: must be less than mu*(1-mu)");
  }

  // Calculate constants a and b based on mu and sigma2
  double cons = ((1 - mu) / sigma2) - (1 / mu);
  double a = cons * (mu * mu);
  double b = cons * mu * (1 - mu);

  // Compute the beta density function for scalar value x
  return R::dbeta(x, a, b, false); // Return beta density for scalar x
}

//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double rbeta_rep_cpp(double mu, double sigma2) {
  double cons = (((1 - mu) / sigma2) - (1 / mu));
  double a = cons * pow(mu, 2);
  double b = cons * mu * (1 - mu);

  // Use R's rbeta function from the Rcpp library
  return R::rbeta(a, b);
}

//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double posteriorvar_cpp(const arma::vec& theta, const arma::mat& adjmatinf, const arma::mat& X, const arma::vec& y,
                        const arma::vec& mutheta, const arma::vec& sigmatheta) {
  double rho = theta[0];
  double psi = theta[1];

  double murho = mutheta[0];
  double mupsi = mutheta[1];

  double sigmarho = sigmatheta[0];
  double sigmapsi = sigmatheta[1];

  // Compute prior densities for rho and psi
  double priorrho = dbeta_rep_cpp(rho, murho, sigmarho);
  double priorpsi = dbeta_rep_cpp(psi, mupsi, sigmapsi);

  // Compute variance-covariance matrix
  arma::mat varcovdagar = vdagar_cpp(adjmatinf, rho, psi);
  arma::mat invvarcovdagar = arma::inv(varcovdagar);

  arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
  arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


  int n = X.n_rows;
  int p = X.n_cols;

  // Compute residuals and variance
  arma::vec residuals = y - X * betaest;
  double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

  // Determinants for posterior calculation
  double detvarcovinv = arma::det(invvarcovdagar);
  double detvarbeta = arma::det(varbeta);

  // Return the posterior variance
  return std::sqrt(detvarcovinv * detvarbeta * std::pow(s2, -(n - p))) * priorrho * priorpsi;
}

//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double posteriorvar2_cpp(const arma::vec& theta, const arma::mat& adjmatinf, const arma::mat& X, const arma::vec& y,
                        const arma::vec& mutheta, const arma::vec& sigmatheta) {
  double rho = theta[0];
  double psi = theta[1];

  double murho = mutheta[0];
  double mupsi = mutheta[1];

  double sigmarho = sigmatheta[0];
  double sigmapsi = sigmatheta[1];

  // Compute prior densities for rho and psi
  double priorrho = dbeta_rep_cpp(rho, murho, sigmarho);
  double priorpsi = dbeta_rep_cpp(psi, mupsi, sigmapsi);

  // Compute variance-covariance matrix
  arma::mat varcovdagar = vdagar_cpp(adjmatinf, rho, psi);
  arma::mat invvarcovdagar = arma::inv(varcovdagar);

  arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
  arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


  int n = X.n_rows;
  int p = X.n_cols;

  // Compute residuals and variance
  arma::vec residuals = y - X * betaest;
  double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

  // Determinants for posterior calculation
  double detvarcovinv = arma::det(invvarcovdagar);
  double detvarbeta = arma::det(varbeta);

  // Return the posterior variance
  return std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho * priorpsi;
}


bool is_sympd(const arma::mat& A) {
  // Check if the matrix is symmetric
  if (!arma::approx_equal(A, arma::trans(A), "absdiff", 1e-8)) {
    return false;
  }

  // Check if the matrix is positive definite by checking if all eigenvalues are positive
  arma::vec eigvals = arma::eig_sym(A);
  return arma::all(eigvals > 0);
}



//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::rowvec mdivide_right_spd_vector(const arma::rowvec& y, const arma::mat& A) {
  // Ensure that A is symmetric and positive definite
  if (A.n_rows != A.n_cols) {
    throw std::invalid_argument("Matrix A must be square.");
  }
  if (!is_sympd(A)) {
    throw std::invalid_argument("Matrix A must be symmetric positive definite.");
  }

  // Cholesky decomposition of A
  arma::mat L = arma::chol(A, "lower");

  // Convert y (row vector) to column vector for the solve process
  arma::colvec y_col = arma::trans(y);  // transpose to column vector

  // Solve the system L * y_col = y (where y_col is the column vector)
  arma::colvec y_solution = arma::solve(L, y_col);

  // Solve the system L^T * x = y_solution
  arma::colvec result = arma::solve(arma::trans(L), y_solution);

  // Return the result as a row vector
  return arma::trans(result);  // convert back to row vector
}


//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat mdivide_right_spd(const arma::mat& A, const arma::mat& B) {
  // Check if A is square and symmetric
  if (A.n_rows != A.n_cols) {
    throw std::invalid_argument("Matrix A must be square.");
  }

  // Perform Cholesky decomposition
  arma::mat L = arma::chol(A); // L is the lower triangular matrix

  // Solve the system L * y = B
  arma::mat y = arma::solve(L, B); // y = L^{-1} * B

  // Now solve L' * x = y for x (which is equivalent to A^{-1} * B)
  arma::mat X = arma::solve(arma::trans(L), y); // X = L'^{-1} * y

  return X; // This returns B * A^{-1}
}

double posterior_spatemp_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                                 const arma::mat& X, const arma::vec& y,
                                 const arma::vec& a, const arma::vec& b);

//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
 // [[Rcpp::export]]
double posterior_spatemp_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                         const arma::mat& X, const arma::vec& y,
                         const arma::vec& a, const arma::vec& b) {
  double rho = theta[0];
  double psi = theta[1];
  double phi = theta[2];

  double arho = a[0];
  double apsi = a[1];
  double aphi = a[2];

  double brho = b[0];
  double bpsi = b[1];
  double bphi = b[2];

  // Compute prior densities for rho and psi
  double priorrho =  R::dbeta(rho, arho, brho, false);
  double priorpsi = R::dbeta(psi, apsi, bpsi,false);
  double priorphi = R::dbeta(phi, aphi, bphi,false);

  // Compute variance-covariance matrix
  mat varcovdagar = spatimecovar_2_cpp(lag,adjmatinf,rho,psi,phi);
  arma::mat invvarcovdagar = arma::inv(varcovdagar);
  arma::mat mult =  X.t() * invvarcovdagar;
  arma::mat varbeta = arma::inv(mult* X);
  arma::vec betaest = varbeta * mult * y;


  int n = X.n_rows;
  int p = X.n_cols;

  // Compute residuals and variance
  arma::vec residuals = y - X * betaest;
  double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

  // Determinants for posterior calculation
  double detvarcovinv = arma::det(invvarcovdagar);
  if (detvarcovinv == R_PosInf) {
    arma::vec detvarcovinvvec(2);
    detvarcovinvvec(0) = 1e300;
    detvarcovinvvec(1) = 1e308;
    detvarcovinv = detvarcovinvvec[R::runif(0, 2)];
  }
  double log_detvarcovinv = log(detvarcovinv);
  double detvarbeta = arma::det(varbeta);
  double log_detvarbeta = log(detvarbeta);
  double s2np = std::pow(s2, -(n - p));
  double log_s2np = -(n-p)*log(s2);

  double posterior_log = 0.5*(log_detvarcovinv + log_detvarbeta +log_s2np) + log(priorrho) +log(priorpsi) +log(priorphi);

  // Return the posterior variance
  //double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho * priorpsi*priorphi;
  return posterior_log;
  //List result;
  //result["detvcovinv"] = detvarcovinv;
  //result["detvarbeta"] = detvarbeta;
  //result["s2expo"] = std::pow(s2, -(n - p));
  //result["posterior"] = posterior;
  //result["priorrho"] = priorrho;
  //result["priorpsi"] = priorpsi;
  //result["priorphi"] = priorphi;
  //return result;
}



double posterior_spatempsar_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                             const arma::mat& X, const arma::vec& y,
                             const arma::vec& a, const arma::vec& b);

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 double posterior_spatempsar_cpp(const arma::vec& theta, const arma::mat& adjmat,int lag,
                              const arma::mat& X, const arma::vec& y,
                              const arma::vec& a, const arma::vec& b) {
   double rho = theta[0];
   double psi = theta[1];
   double phi = theta[2];

   double arho = a[0];
   double apsi = a[1];
   double aphi = a[2];

   double brho = b[0];
   double bpsi = b[1];
   double bphi = b[2];

   // Compute prior densities for rho and psi
   double priorrho =  R::dbeta(rho, arho, brho, false);
   double priorpsi = R::dbeta(psi, apsi, bpsi,false);
   double priorphi = R::dbeta(phi, aphi, bphi,false);

   // Compute variance-covariance matrix
   mat varcovdagar = spatimecovarsar_2_cpp(lag,adjmat,rho,psi,phi);
   arma::mat invvarcovdagar = arma::inv(varcovdagar);
   arma::mat mult =  X.t() * invvarcovdagar;
   arma::mat varbeta = arma::inv(mult* X);
   arma::vec betaest = varbeta * mult * y;


   int n = X.n_rows;
   int p = X.n_cols;

   // Compute residuals and variance
   arma::vec residuals = y - X * betaest;
   double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

   // Determinants for posterior calculation
   double detvarcovinv = arma::det(invvarcovdagar);
   if (detvarcovinv == R_PosInf) {
     arma::vec detvarcovinvvec(2);
     detvarcovinvvec(0) = 1e300;
     detvarcovinvvec(1) = 1e308;
     detvarcovinv = detvarcovinvvec[R::runif(0, 2)];
   }
   double log_detvarcovinv = log(detvarcovinv);
   double detvarbeta = arma::det(varbeta);
   double log_detvarbeta = log(detvarbeta);
   double s2np = std::pow(s2, -(n - p));
   double log_s2np = -(n-p)*log(s2);

   double posterior_log = 0.5*(log_detvarcovinv + log_detvarbeta +log_s2np) + log(priorrho) +log(priorpsi) +log(priorphi);

   // Return the posterior variance
   //double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho * priorpsi*priorphi;
   return posterior_log;
   //List result;
   //result["detvcovinv"] = detvarcovinv;
   //result["detvarbeta"] = detvarbeta;
   //result["s2expo"] = std::pow(s2, -(n - p));
   //result["posterior"] = posterior;
   //result["priorrho"] = priorrho;
   //result["priorpsi"] = priorpsi;
   //result["priorphi"] = priorphi;
   //return result;
 }




// [[Rcpp::export]]
double dbeta_ab(double y, double alpha, double beta,
                double a, double b, bool log_out = false) {
  if (y <= a || y >= b) {
    return log_out ? R_NegInf : 0.0;
  }

  // scale to (0,1)
  double z = (y - a) / (b - a);

  // base beta density
  double val = R::dbeta(z, alpha, beta, log_out);

  // adjust for transformation
  if (log_out) {
    val -= std::log(b - a);
  } else {
    val /= (b - a);
  }

  return val;
}


Rcpp::List posterior_spatempcens_2_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                                     const arma::mat& X, const arma::vec& y,
                                     const arma::vec& a, const arma::vec& b);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 Rcpp::List posterior_spatempcens_2_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                                      const arma::mat& X, const arma::vec& y,
                                      const arma::vec& a, const arma::vec& b) {
   double rho = theta[0];
   double psi = theta[1];
   double phi = theta[2];

   double arho = a[0];
   double apsi = a[1];
   double aphi = a[2];

   double brho = b[0];
   double bpsi = b[1];
   double bphi = b[2];

   // Compute prior densities for rho and psi
   double priorrho =  R::dbeta(rho, arho, brho, false);
   double priorpsi = R::dbeta(psi, apsi, bpsi,false);
   double priorphi = dbeta_ab(phi,aphi, bphi,-1.0, 1.0, false);

   // Compute variance-covariance matrix
   mat varcovdagar = spatimecovar_2_cpp(lag,adjmatinf,rho,psi,phi);
   arma::mat invvarcovdagar = arma::inv(varcovdagar);

   arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
   arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


   int n = X.n_rows;
   int p = X.n_cols;

   // Compute residuals and variance
   arma::vec residuals = y - X * betaest;
   double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

   // Determinants for posterior calculation
   double detvarcovinv = arma::det(invvarcovdagar);
   if (detvarcovinv == R_PosInf) {
     detvarcovinv = 1e308;
   }
   double log_detvarcovinv = log(detvarcovinv);
   double detvarbeta = arma::det(varbeta);
   double log_detvarbeta = log(detvarbeta);
   double s2np = std::pow(s2, -(n - p));
   double log_s2np = -(n-p)*log(s2);

   // Return the posterior variance
   double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho * priorpsi*priorphi;
   double posterior_log = 0.5*(log_detvarcovinv + log_detvarbeta +log_s2np) + log(priorrho) +log(priorpsi) +log(priorphi);
   return Rcpp::List::create(
     Rcpp::Named("priorrho") = priorrho,
     Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("priorphi") = priorphi,
     Rcpp::Named("varcovspatemp") = varcovdagar,
     Rcpp::Named("varbeta") = varbeta,
     Rcpp::Named("betaest") = betaest,
     Rcpp::Named("S2") = s2,
     //Rcpp::Named("priorphi") = priorphi,
     //Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("posterior") = posterior,
     Rcpp::Named("posterior_log") = posterior_log
   );
   //List result;
   //result["detvcovinv"] = detvarcovinv;
   //result["detvarbeta"] = detvarbeta;
   //result["s2expo"] = std::pow(s2, -(n - p));
   //result["posterior"] = posterior;
   //result["priorrho"] = priorrho;
   //result["priorpsi"] = priorpsi;
   //result["priorphi"] = priorphi;
   //return result;
 }


Rcpp::List posterior_spatempcenssar_2_cpp(const arma::vec& theta, const arma::mat& adjmat,int lag,
                                       const arma::mat& X, const arma::vec& y,
                                       const arma::vec& a, const arma::vec& b);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 Rcpp::List posterior_spatempcenssar_2_cpp(const arma::vec& theta, const arma::mat& adjmat,int lag,
                                        const arma::mat& X, const arma::vec& y,
                                        const arma::vec& a, const arma::vec& b) {
   double rho = theta[0];
   double psi = theta[1];
   double phi = theta[2];

   double arho = a[0];
   double apsi = a[1];
   double aphi = a[2];

   double brho = b[0];
   double bpsi = b[1];
   double bphi = b[2];

   // Compute prior densities for rho and psi
   double priorrho =  dbeta_ab(rho, arho, brho, -1.0, 1.0, false);
   double priorpsi = R::dbeta(psi, apsi, bpsi,false);
   double priorphi = dbeta_ab(phi,aphi, bphi,-1.0, 1.0, false);

   // Compute variance-covariance matrix
   mat varcovdagar = spatimecovarsar_2_cpp(lag,adjmat,rho,psi,phi);
   arma::mat invvarcovdagar = arma::inv(varcovdagar);

   arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
   arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


   int n = X.n_rows;
   int p = X.n_cols;

   // Compute residuals and variance
   arma::vec residuals = y - X * betaest;
   double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

   // Determinants for posterior calculation
   double detvarcovinv = arma::det(invvarcovdagar);
   if (detvarcovinv == R_PosInf) {
     detvarcovinv = 1e308;
   }
   double log_detvarcovinv = log(detvarcovinv);
   double detvarbeta = arma::det(varbeta);
   double log_detvarbeta = log(detvarbeta);
   double s2np = std::pow(s2, -(n - p));
   double log_s2np = -(n-p)*log(s2);

   // Return the posterior variance
   double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho * priorpsi*priorphi;
   double posterior_log = 0.5*(log_detvarcovinv + log_detvarbeta +log_s2np) + log(priorrho) +log(priorpsi) +log(priorphi);
   return Rcpp::List::create(
     Rcpp::Named("priorrho") = priorrho,
     Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("priorphi") = priorphi,
     Rcpp::Named("varcovspatemp") = varcovdagar,
     Rcpp::Named("varbeta") = varbeta,
     Rcpp::Named("betaest") = betaest,
     Rcpp::Named("S2") = s2,
     //Rcpp::Named("priorphi") = priorphi,
     //Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("posterior") = posterior,
     Rcpp::Named("posterior_log") = posterior_log
   );
   ////List result;
   //result["detvcovinv"] = detvarcovinv;
   //result["detvarbeta"] = detvarbeta;
   //result["s2expo"] = std::pow(s2, -(n - p));
   //result["posterior"] = posterior;
   //result["priorrho"] = priorrho;
   //result["priorpsi"] = priorpsi;
   ////result["priorphi"] = priorphi;
   ////return result;
 }







Rcpp::List posterior_spatempcens_ar2_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                                       const arma::mat& X, const arma::vec& y,
                                       const arma::vec& a, const arma::vec& b);

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 Rcpp::List posterior_spatempcens_ar2_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                                        const arma::mat& X, const arma::vec& y,
                                        const arma::vec& a, const arma::vec& b) {
   double rho = theta[0];
   double psi = theta[1];
   double phi1 = theta[2];
   double phi2 = theta[3];

   double arho = a[0];
   double apsi = a[1];
   double aphi1 = a[2];
   double aphi2 = a[3];

   double brho = b[0];
   double bpsi = b[1];
   double bphi1 = b[2];
   double bphi2 = b[3];

   // Compute prior densities for rho and psi
   double priorrho =  R::dbeta(rho, arho, brho, false);
   double priorpsi = R::dbeta(psi, apsi, bpsi,false);
   double priorphi2 = dbeta_ab(phi2,aphi2, bphi2,-1.0, 1.0, false);
   double priorphi1 = dbeta_ab(phi1,aphi1, bphi1, phi2 - 1.0, 1.0-phi2, false);

   // Compute variance-covariance matrix
   mat varcovdagar = spatimecovar_ar2_cpp(lag,adjmatinf,rho,psi,phi1,phi2);
   arma::mat invvarcovdagar = arma::inv(varcovdagar);

   arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
   arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


   int n = X.n_rows;
   int p = X.n_cols;

   // Compute residuals and variance
   arma::vec residuals = y - X * betaest;
   double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

   // Determinants for posterior calculation
   double detvarcovinv = arma::det(invvarcovdagar);
   if (detvarcovinv == R_PosInf) {
     detvarcovinv = 1e308;
   }
   double log_detvarcovinv = log(detvarcovinv);
   double detvarbeta = arma::det(varbeta);
   double log_detvarbeta = log(detvarbeta);
   double s2np = std::pow(s2, -(n - p));
   double log_s2np = -(n-p)*log(s2);

   // Return the posterior variance
   double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho * priorpsi*priorphi1*priorphi2;
   double posterior_log = 0.5*(log_detvarcovinv + log_detvarbeta +log_s2np) + log(priorrho) +log(priorpsi) +log(priorphi1) +log(priorphi2);
   return Rcpp::List::create(
     Rcpp::Named("priorrho") = priorrho,
     Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("priorphi1") = priorphi1,
     Rcpp::Named("priorphi1") = priorphi2,
     Rcpp::Named("varcovspatemp") = varcovdagar,
     Rcpp::Named("varbeta") = varbeta,
     Rcpp::Named("betaest") = betaest,
     Rcpp::Named("S2") = s2,
     //Rcpp::Named("priorphi") = priorphi,
     //Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("posterior") = posterior,
     Rcpp::Named("posterior_log") = posterior_log
   );
   //List result;
   //result["detvcovinv"] = detvarcovinv;
   //result["detvarbeta"] = detvarbeta;
   //result["s2expo"] = std::pow(s2, -(n - p));
   //result["posterior"] = posterior;
   //result["priorrho"] = priorrho;
   //result["priorpsi"] = priorpsi;
   //result["priorphi"] = priorphi;
   //return result;
 }


Rcpp::List posterior_spatempcenssar_ar2_cpp(const arma::vec& theta, const arma::mat& adjmat,int lag,
                                         const arma::mat& X, const arma::vec& y,
                                         const arma::vec& a, const arma::vec& b);

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 Rcpp::List posterior_spatempcenssar_ar2_cpp(const arma::vec& theta, const arma::mat& adjmat,int lag,
                                          const arma::mat& X, const arma::vec& y,
                                          const arma::vec& a, const arma::vec& b) {
   double rho = theta[0];
   double psi = theta[1];
   double phi1 = theta[2];
   double phi2 = theta[3];

   double arho = a[0];
   double apsi = a[1];
   double aphi1 = a[2];
   double aphi2 = a[3];

   double brho = b[0];
   double bpsi = b[1];
   double bphi1 = b[2];
   double bphi2 = b[3];

   // Compute prior densities for rho and psi
   double priorrho =  dbeta_ab(rho, arho, brho, -1.0,1.0, false);
   double priorpsi = R::dbeta(psi, apsi, bpsi,false);
   double priorphi2 = dbeta_ab(phi2,aphi2, bphi2,-1.0, 1.0, false);
   double priorphi1 = dbeta_ab(phi1,aphi1, bphi1, phi2 - 1.0, 1.0-phi2, false);

   // Compute variance-covariance matrix
   mat varcovdagar = spatimecovarsar_ar2_cpp(lag,adjmat,rho,psi,phi1,phi2);
   arma::mat invvarcovdagar = arma::inv(varcovdagar);

   arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
   arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


   int n = X.n_rows;
   int p = X.n_cols;

   // Compute residuals and variance
   arma::vec residuals = y - X * betaest;
   double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

   // Determinants for posterior calculation
   double detvarcovinv = arma::det(invvarcovdagar);
   if (detvarcovinv == R_PosInf) {
     detvarcovinv = 1e308;
   }
   double log_detvarcovinv = log(detvarcovinv);
   double detvarbeta = arma::det(varbeta);
   double log_detvarbeta = log(detvarbeta);
   double s2np = std::pow(s2, -(n - p));
   double log_s2np = -(n-p)*log(s2);

   // Return the posterior variance
   double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho * priorpsi*priorphi1*priorphi2;
   double posterior_log = 0.5*(log_detvarcovinv + log_detvarbeta +log_s2np) + log(priorrho) +log(priorpsi) +log(priorphi1) +log(priorphi2);
   return Rcpp::List::create(
     Rcpp::Named("priorrho") = priorrho,
     Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("priorphi1") = priorphi1,
     Rcpp::Named("priorphi1") = priorphi2,
     Rcpp::Named("varcovspatemp") = varcovdagar,
     Rcpp::Named("varbeta") = varbeta,
     Rcpp::Named("betaest") = betaest,
     Rcpp::Named("S2") = s2,
     //Rcpp::Named("priorphi") = priorphi,
     //Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("posterior") = posterior,
     Rcpp::Named("posterior_log") = posterior_log
   );
   //List result;
   //result["detvcovinv"] = detvarcovinv;
   //result["detvarbeta"] = detvarbeta;
   //result["s2expo"] = std::pow(s2, -(n - p));
   //result["posterior"] = posterior;
   //result["priorrho"] = priorrho;
   //result["priorpsi"] = priorpsi;
   //result["priorphi"] = priorphi;
   //return result;
 }


Rcpp::List posterior_spatempcens_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                             const arma::mat& X, const arma::vec& y,
                             const arma::vec& a, const arma::vec& b);
//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
Rcpp::List posterior_spatempcens_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                              const arma::mat& X, const arma::vec& y,
                              const arma::vec& a, const arma::vec& b) {
   double rho = theta[0];
   double psi = theta[1];
   double phi = theta[2];

   double arho = a[0];
   double apsi = a[1];
   double aphi = a[2];

   double brho = b[0];
   double bpsi = b[1];
   double bphi = b[2];

   // Compute prior densities for rho and psi
   double priorrho =  R::dbeta(rho, arho, brho, false);
   double priorpsi = R::dbeta(psi, apsi, bpsi,false);
   double priorphi = R::dbeta(phi, aphi, bphi,false);

   // Compute variance-covariance matrix
   mat varcovdagar = spatimecovar_2_cpp(lag,adjmatinf,rho,psi,phi);
   arma::mat invvarcovdagar = arma::inv(varcovdagar);

   arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
   arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


   int n = X.n_rows;
   int p = X.n_cols;

   // Compute residuals and variance
   arma::vec residuals = y - X * betaest;
   double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

   // Determinants for posterior calculation
   double detvarcovinv = arma::det(invvarcovdagar);
   if (detvarcovinv == R_PosInf) {
     detvarcovinv = 1e308;
   }
   double log_detvarcovinv = log(detvarcovinv);
   double detvarbeta = arma::det(varbeta);
   double log_detvarbeta = log(detvarbeta);
   double s2np = std::pow(s2, -(n - p));
   double log_s2np = -(n-p)*log(s2);

   // Return the posterior variance
   double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho * priorpsi*priorphi;
   double posterior_log = 0.5*(log_detvarcovinv + log_detvarbeta +log_s2np) + log(priorrho) +log(priorpsi) +log(priorphi);
   return Rcpp::List::create(
     Rcpp::Named("priorrho") = priorrho,
     Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("priorphi") = priorphi,
     Rcpp::Named("varcovspatemp") = varcovdagar,
     Rcpp::Named("varbeta") = varbeta,
     Rcpp::Named("betaest") = betaest,
     Rcpp::Named("S2") = s2,
     //Rcpp::Named("priorphi") = priorphi,
     //Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("posterior") = posterior,
     Rcpp::Named("posterior_log") = posterior_log
   );
   //List result;
   //result["detvcovinv"] = detvarcovinv;
   //result["detvarbeta"] = detvarbeta;
   //result["s2expo"] = std::pow(s2, -(n - p));
   //result["posterior"] = posterior;
   //result["priorrho"] = priorrho;
   //result["priorpsi"] = priorpsi;
   //result["priorphi"] = priorphi;
   //return result;
 }




Rcpp::List posterior_spatempcenscar_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                                     const arma::mat& X, const arma::vec& y,
                                     const arma::vec& a, const arma::vec& b);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 Rcpp::List posterior_spatempcenscar_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                                      const arma::mat& X, const arma::vec& y,
                                      const arma::vec& a, const arma::vec& b) {
   double rho = theta[0];
   double psi = theta[1];
   double phi = theta[2];

   double arho = a[0];
   double apsi = a[1];
   double aphi = a[2];

   double brho = b[0];
   double bpsi = b[1];
   double bphi = b[2];

   // Compute prior densities for rho and psi
   double priorrho =  R::dbeta(rho, arho, brho, false);
   double priorpsi = R::dbeta(psi, apsi, bpsi,false);
   double priorphi = R::dbeta(phi, aphi, bphi,false);

   // Compute variance-covariance matrix
   mat varcovdagar = spatimecovarcar_2_cpp(lag,adjmatinf,rho,psi,phi);
   arma::mat invvarcovdagar = arma::inv(varcovdagar);

   arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
   arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


   int n = X.n_rows;
   int p = X.n_cols;

   // Compute residuals and variance
   arma::vec residuals = y - X * betaest;
   double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

   // Determinants for posterior calculation
   double detvarcovinv = arma::det(invvarcovdagar);
   if (detvarcovinv == R_PosInf) {
     detvarcovinv = 1e308;
   }
   double log_detvarcovinv = log(detvarcovinv);
   double detvarbeta = arma::det(varbeta);
   double log_detvarbeta = log(detvarbeta);
   double s2np = std::pow(s2, -(n - p));
   double log_s2np = -(n-p)*log(s2);

   // Return the posterior variance
   double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho * priorpsi*priorphi;
   double posterior_log = 0.5*(log_detvarcovinv + log_detvarbeta +log_s2np) + log(priorrho) +log(priorpsi) +log(priorphi);
   return Rcpp::List::create(
     Rcpp::Named("priorrho") = priorrho,
     Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("priorphi") = priorphi,
     Rcpp::Named("varcovspatemp") = varcovdagar,
     Rcpp::Named("varbeta") = varbeta,
     Rcpp::Named("betaest") = betaest,
     Rcpp::Named("S2") = s2,
     //Rcpp::Named("priorphi") = priorphi,
     //Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("posterior") = posterior,
     Rcpp::Named("posterior_log") = posterior_log
   );
   //List result;
   //result["detvcovinv"] = detvarcovinv;
   //result["detvarbeta"] = detvarbeta;
   //result["s2expo"] = std::pow(s2, -(n - p));
   //result["posterior"] = posterior;
   //result["priorrho"] = priorrho;
   //result["priorpsi"] = priorpsi;
   //result["priorphi"] = priorphi;
   //return result;
 }



















//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
// [[Rcpp::export]]
arma::mat vector_to_matrix(const arma::vec& v, int nrows, int ncols) {
  // Ensure size matches
  if (v.n_elem != nrows * ncols) {
    stop("Vector length does not match specified matrix dimensions.");
  }
  // Reshape (Column-major order)
  return arma::mat(v).reshape(nrows, ncols);
}




Rcpp::List posterior_spatempcens_zero_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                                    const arma::mat& X, const arma::vec& y,
                                     const arma::vec& a, const arma::vec& b);
//'@useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 Rcpp::List posterior_spatempcens_zero_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                                      const arma::mat& X, const arma::vec& y,
                                      const arma::vec& a, const arma::vec& b) {
   double rho = theta[0];
   double phi = theta[1];

   double arho = a[0];
   double aphi = a[1];

   double brho = b[0];
   double bphi = b[1];

   // Compute prior densities for rho and psi
   double priorrho =  R::dbeta(rho, arho, brho, false);
   double priorphi = dbeta_ab(phi, aphi, bphi,-1.0,1.0, false);
   // Compute variance-covariance matrix
   mat varcovspadagar = varcovspadagar_cpp(adjmatinf,rho);
   mat varcovtime = ar1_correlation_matrix_2(phi,lag);
   mat varcovdagar = arma::kron(varcovspadagar,varcovtime);
   arma::mat invvarcovspadagar = arma::inv(varcovspadagar);
   arma::mat invvarcovtime = arma::inv(varcovtime);
   arma::mat invvarcovdagar = arma::kron(invvarcovspadagar,invvarcovtime);
   arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
   arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


   int n = X.n_rows;
   int p = X.n_cols;
   int nspa = varcovspadagar.n_rows;
   int ntime = varcovtime.n_rows;

   // Compute residuals and variance
   arma::vec residuals = y - X * betaest;
   arma::mat residualsmat = vector_to_matrix(residuals,ntime,nspa);
   arma::mat matquad = invvarcovspadagar*residualsmat.t()*invvarcovtime*residualsmat;
   double s2 = arma::trace(matquad)/(n-p);
   //double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

   // Determinants for posterior calculation
   double detvarcovspainv = arma::det(invvarcovdagar);
   double detvarcovtimeinv = arma::det(invvarcovtime);
   double detvarcovinv = std::pow(detvarcovspainv,nspa)*std::pow(detvarcovtimeinv,ntime);
   if (detvarcovspainv == R_PosInf) {
     detvarcovspainv = 1e308;
   }
   if (detvarcovtimeinv == R_PosInf) {
     detvarcovtimeinv = 1e308;
   }

   if (detvarcovinv == R_PosInf) {
     detvarcovinv = 1e308;
   }

   double log_detvarcovinv = log(detvarcovinv);
   double detvarbeta = arma::det(varbeta);
   double log_detvarbeta = log(detvarbeta);
   double s2np = std::pow(s2, -(n - p));
   double log_s2np = -(n-p)*log(s2);

   // Return the posterior variance
   double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho *priorphi;
   double posterior_log = 0.5*(log_detvarcovinv + log_detvarbeta +log_s2np) + log(priorrho) +log(priorphi);
   return Rcpp::List::create(
     Rcpp::Named("priorrho") = priorrho,
     Rcpp::Named("priorphi") = priorphi,
     Rcpp::Named("varcovspatemp") = varcovdagar,
     Rcpp::Named("varbeta") = varbeta,
     Rcpp::Named("betaest") = betaest,
     Rcpp::Named("S2") = s2,
     //Rcpp::Named("priorphi") = priorphi,
     //Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("posterior") = posterior,
     Rcpp::Named("posterior_log") = posterior_log
   );
   //List result;
   //result["detvcovinv"] = detvarcovinv;
   //result["detvarbeta"] = detvarbeta;
   //result["s2expo"] = std::pow(s2, -(n - p));
   //result["posterior"] = posterior;
   //result["priorrho"] = priorrho;
   //result["priorpsi"] = priorpsi;
   //result["priorphi"] = priorphi;
   //return result;
 }



Rcpp::List posterior_spatempcens_ar2_zero_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                                          const arma::mat& X, const arma::vec& y,
                                          const arma::vec& a, const arma::vec& b);
//'@useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 Rcpp::List posterior_spatempcens_ar2_zero_cpp(const arma::vec& theta, const arma::mat& adjmatinf,int lag,
                                           const arma::mat& X, const arma::vec& y,
                                           const arma::vec& a, const arma::vec& b) {
   double rho = theta[0];
   double phi1 = theta[1];
   double phi2 = theta[2];

   double arho = a[0];
   double aphi1 = a[1];
   double aphi2 = a[1];

   double brho = b[0];
   double bphi1 = b[1];
   double bphi2 = b[2];

   // Compute prior densities for rho and psi
   double priorrho =  R::dbeta(rho, arho, brho, false);
   double priorphi2 = dbeta_ab(phi2,aphi2, bphi2,-1.0, 1.0, false);
   double priorphi1 = dbeta_ab(phi1,aphi1, bphi1, phi2-1.0, 1.0-phi2, false);

   // Compute variance-covariance matrix
   mat varcovspadagar = varcovspadagar_cpp(adjmatinf,rho);
   mat varcovtime = corr_ar2_fast(lag,phi1,phi2,true);
   mat varcovdagar = arma::kron(varcovspadagar,varcovtime);
   arma::mat invvarcovspadagar = arma::inv(varcovspadagar);
   arma::mat invvarcovtime = arma::inv(varcovtime);
   arma::mat invvarcovdagar = arma::kron(invvarcovspadagar,invvarcovtime);
   arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
   arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


   int n = X.n_rows;
   int p = X.n_cols;
   int nspa = varcovspadagar.n_rows;
   int ntime = varcovtime.n_rows;

   // Compute residuals and variance
   arma::vec residuals = y - X * betaest;
   arma::mat residualsmat = vector_to_matrix(residuals,ntime,nspa);
   arma::mat matquad = invvarcovspadagar*residualsmat.t()*invvarcovtime*residualsmat;
   double s2 = arma::trace(matquad)/(n-p);
   //double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

   // Determinants for posterior calculation
   double detvarcovspainv = arma::det(invvarcovdagar);
   double detvarcovtimeinv = arma::det(invvarcovtime);
   double detvarcovinv = std::pow(detvarcovspainv,nspa)*std::pow(detvarcovtimeinv,ntime);
   if (detvarcovspainv == R_PosInf) {
     detvarcovspainv = 1e308;
   }
   if (detvarcovtimeinv == R_PosInf) {
     detvarcovtimeinv = 1e308;
   }

   if (detvarcovinv == R_PosInf) {
     detvarcovinv = 1e308;
   }

   double log_detvarcovinv = log(detvarcovinv);
   double detvarbeta = arma::det(varbeta);
   double log_detvarbeta = log(detvarbeta);
   double s2np = std::pow(s2, -(n - p));
   double log_s2np = -(n-p)*log(s2);

   // Return the posterior variance
   double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho *priorphi1*priorphi2;
   double posterior_log = 0.5*(log_detvarcovinv + log_detvarbeta +log_s2np) + log(priorrho) +log(priorphi1) + log(priorphi2);
   return Rcpp::List::create(
     Rcpp::Named("priorrho") = priorrho,
     Rcpp::Named("priorphi1") = priorphi1,
     Rcpp::Named("priorphi2") = priorphi2,
     Rcpp::Named("varcovspatemp") = varcovdagar,
     Rcpp::Named("varbeta") = varbeta,
     Rcpp::Named("betaest") = betaest,
     Rcpp::Named("S2") = s2,
     //Rcpp::Named("priorphi") = priorphi,
     //Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("posterior") = posterior,
     Rcpp::Named("posterior_log") = posterior_log
   );
   //List result;
   //result["detvcovinv"] = detvarcovinv;
   //result["detvarbeta"] = detvarbeta;
   //result["s2expo"] = std::pow(s2, -(n - p));
   //result["posterior"] = posterior;
   //result["priorrho"] = priorrho;
   //result["priorpsi"] = priorpsi;
   //result["priorphi"] = priorphi;
   //return result;
 }



//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double posteriorvar3_cpp(const arma::vec& theta, const arma::mat& adjmatinf, const arma::mat& X, const arma::vec& y,
                         const arma::vec& a, const arma::vec& b) {
  double rho = theta[0];
  double psi = theta[1];

  double arho = a[0];
  double apsi = a[1];

  double brho = b[0];
  double bpsi = b[1];

  // Compute prior densities for rho and psi
  double priorrho =  R::dbeta(rho, arho, brho, false);
  double priorpsi = R::dbeta(psi, apsi, bpsi,false);

  // Compute variance-covariance matrix
  arma::mat varcovdagar = vdagar_cpp(adjmatinf, rho, psi);
  arma::mat invvarcovdagar = arma::inv(varcovdagar);

  arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
  arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


  int n = X.n_rows;
  int p = X.n_cols;

  // Compute residuals and variance
  arma::vec residuals = y - X * betaest;
  double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

  // Determinants for posterior calculation
  double detvarcovinv = arma::det(invvarcovdagar);
  double detvarbeta = arma::det(varbeta);

  // Return the posterior variance
  return std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho * priorpsi;
}
// [[Rcpp::export]]
Rcpp::List posteriorvar4_cpp(const arma::vec& theta, const arma::mat& adjmatinf, const arma::mat& X, const arma::vec& y,
                         const arma::vec& a, const arma::vec& b) {
  double rho = theta[0];
  double psi = theta[1];

  double arho = a[0];
  double apsi = a[1];

  double brho = b[0];
  double bpsi = b[1];

  // Compute prior densities for rho and psi
  double priorrho =  R::dbeta(rho, arho, brho, false);
  double priorpsi = R::dbeta(psi, apsi, bpsi,false);

  // Compute variance-covariance matrix
  arma::mat varcovdagar = vdagar_cpp(adjmatinf, rho, psi);
  arma::mat invvarcovdagar = arma::inv(varcovdagar);

  arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
  arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


  int n = X.n_rows;
  int p = X.n_cols;

  // Compute residuals and variance
  arma::vec residuals = y - X * betaest;
  double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

  // Determinants for posterior calculation
  double detvarcovinv = arma::det(invvarcovdagar);
  double detvarbeta = arma::det(varbeta);

  // Return the posterior variance
  double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho * priorpsi;

  return Rcpp::List::create(
    Rcpp::Named("posterior") = posterior,
    Rcpp::Named("invvarcov") = invvarcovdagar
  );
}











//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
Rcpp::List compute_exponentiated_matrices(const arma::mat& A, int n) {
  int size = A.n_rows;
  Rcpp::List matrices(n);  // Create a list to store the results

  // Initialize the first matrix (A^1)
  arma::mat current = A;
  matrices[0] = current;  // Store A^1

  // Compute matrices A^2, A^3, ..., A^(n-1)
  for (int i = 1; i < n; ++i) {
    current = current * A;  // Compute A^i
    matrices[i] = current;  // Store A^i
  }

  return matrices;  // Return the list of matrices
}

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
Rcpp::List compute_minima_and_argmin(const Rcpp::List& matrices) {
  if (matrices.size() == 0) {
    throw std::invalid_argument("The list of matrices is empty.");
  }

  // Initialize the minima matrix with the first matrix and create an index matrix
  arma::mat min_matrix = Rcpp::as<arma::mat>(matrices[0]);
  arma::mat index_matrix(min_matrix.n_rows, min_matrix.n_cols, arma::fill::zeros);

  // Set all elements of the index matrix to 1 (corresponding to the first matrix, A^1)
  index_matrix.fill(1);

  // Loop through all matrices to compute the element-wise minima and track the indices
  for (size_t i = 1; i < static_cast<size_t>(matrices.size()); ++i) {  // Start from 1, as we've already set the first matrix
    arma::mat current_matrix = Rcpp::as<arma::mat>(matrices[i]);  // Convert list element to matrix

    // Find the positions where the current matrix is less than the min_matrix
    arma::umat less_than_mask = (current_matrix < min_matrix);

    // Extract the positions where the current matrix is less than the min_matrix
    arma::uvec positions = arma::find(less_than_mask);

    // Update min_matrix with new minimum values where the mask is true
    min_matrix.elem(positions) = current_matrix.elem(positions);

    // Create a temporary vector filled with the current exponent (i + 1)
    arma::vec index_update(positions.n_elem);  // Now positions.n_elem works as positions is a uvec
    index_update.fill(i + 1);

    // Update index_matrix with the current exponent where the minimum occurred
    index_matrix.elem(positions) = index_update;
  }

  // Create an R list to return the results
  return Rcpp::List::create(
    Rcpp::Named("minima_matrix") = min_matrix,
    Rcpp::Named("index_matrix") = index_matrix
  );
}

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat compute_argmin_only(const Rcpp::List& matrices) {
  if (matrices.size() == 0) {
    throw std::invalid_argument("The list of matrices is empty.");
  }

  // Initialize the minima matrix with the first matrix and create an index matrix
  arma::mat min_matrix = Rcpp::as<arma::mat>(matrices[0]);
  arma::mat index_matrix(min_matrix.n_rows, min_matrix.n_cols, arma::fill::zeros);

  // Set all elements of the index matrix to 1 (corresponding to the first matrix, A^1)
  index_matrix.fill(1);

  // Loop through all matrices to compute the element-wise minima and track the indices
  for (size_t i = 1; i < static_cast<size_t>(matrices.size()); ++i) {  // Start from 1, as we've already set the first matrix
    arma::mat current_matrix = Rcpp::as<arma::mat>(matrices[i]);  // Convert list element to matrix

    // Find the positions where the current matrix is less than the min_matrix
    arma::umat less_than_mask = (current_matrix < min_matrix);

    // Extract the positions where the current matrix is less than the min_matrix
    arma::uvec positions = arma::find(less_than_mask);

    // Update the index_matrix with the current exponent where the minimum occurred
    arma::vec index_update(positions.n_elem);
    index_update.fill(i + 1);
    index_matrix.elem(positions) = index_update;

    // We don't need to update the min_matrix, as we are only interested in the index matrix
  }

  // Return the index matrix
  return index_matrix;
}

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
Rcpp::List compute_squared_exponentiated_matrices(const arma::mat& A, int n) {
  int size = A.n_rows;
  Rcpp::List squared_matrices(n);  // Preallocate list

  arma::mat current = A;
  arma::mat squared = current * current;

  squared_matrices[0] = squared;  // Store (A^1)^2

  // Use matrix exponentiation and squaring more efficiently
  for (int i = 2; i <= n; ++i) {
    current *= A;                // Compute A^i (matrix exponentiation)
    squared = current * current;  // Compute (A^i)^2
    squared_matrices[i-1] = squared;  // Store (A^i)^2
  }

  return squared_matrices;
}

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat compute_argmin_only_2(const Rcpp::List& matrices) {
  int n = matrices.size();
  arma::mat first_matrix = Rcpp::as<arma::mat>(matrices[0]);
  int rows = first_matrix.n_rows;
  int cols = first_matrix.n_cols;

  // Preallocate matrices for argmin indices and minimum values
  arma::mat argmin_matrix(rows, cols, arma::fill::ones); // Start with index 1
  arma::mat min_values = first_matrix; // Initialize with the first matrix values

  // Efficiently find the minimum index for each element
  for (int k = 1; k < n; ++k) {
    arma::mat next_matrix = Rcpp::as<arma::mat>(matrices[k]);

    // Update the argmin matrix using Armadillo element-wise comparisons
    argmin_matrix.elem(arma::find(next_matrix < min_values)).fill(k + 1);  // Store k+1 as 1-based index
    min_values = arma::min(min_values, next_matrix);  // Update the minimum values matrix
  }

  return argmin_matrix;
}

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double difference(int i, int j, const Rcpp::List& WFval, int n) {

  // Adjust i and j for zero-based indexing
  int idx_i = i - 1;
  int idx_j = j - 1;

  // Ensure the indices are valid
  if (idx_i < 0 || idx_i >= WFval.size() || idx_j < 0 || idx_j >= WFval.size()) {
    Rcpp::stop("Index out of bounds");
  }

  // Extract matrices and normalize them row-wise
  arma::mat WF1star = Rcpp::as<arma::mat>(WFval[idx_i]);
  arma::mat WF2star = Rcpp::as<arma::mat>(WFval[idx_j]);

  // Normalize each row: sum each row and transpose the sum to match dimensions
  WF1star.each_row() /= arma::sum(WF1star, 1).t();  // Correct the dimension by transposing
  WF2star.each_row() /= arma::sum(WF2star, 1).t();  // Correct the dimension by transposing

  // Compute squared exponentiated matrices up to order n-1
  Rcpp::List p1 = compute_squared_exponentiated_matrices(WF1star, n-1);
  Rcpp::List p2 = compute_squared_exponentiated_matrices(WF2star, n-1);

  // Compute argmin matrices for both WF1star and WF2star
  arma::mat path1 = compute_argmin_only_2(p1);
  arma::mat path2 = compute_argmin_only_2(p2);

  // Compute squared difference between path matrices and sum up
  double diff = arma::accu(arma::pow(path1 - path2, 2));

  return diff;
}

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat generate_random_spanning_tree_optimized(const arma::mat& adj_matrix) {
  int n = adj_matrix.n_rows;  // Number of nodes
  arma::mat tree_adj_matrix = arma::zeros<arma::mat>(n, n);  // Initialize spanning tree adjacency matrix

  // Random number generation
  std::random_device rd;
  std::mt19937 gen(rd());

  // Tracking visited nodes
  std::vector<bool> visited(n, false);
  std::vector<int> stack;  // Stack for DFS

  // Start at a random node
  std::uniform_int_distribution<> distrib(0, n-1);
  int start_node = distrib(gen);
  visited[start_node] = true;
  stack.push_back(start_node);

  int edges_added = 0;

  // While the tree has less than n-1 edges
  while (edges_added < n - 1) {
    if (stack.empty()) break;  // Safety check (should not happen for connected graphs)

    int current_node = stack.back();
    stack.pop_back();

    // Gather all unvisited neighbors
    std::vector<int> neighbors;
    for (int j = 0; j < n; ++j) {
      if (adj_matrix(current_node, j) == 1 && !visited[j]) {
        neighbors.push_back(j);
      }
    }

    if (!neighbors.empty()) {
      // Shuffle neighbors to ensure randomness
      std::shuffle(neighbors.begin(), neighbors.end(), gen);

      // Select the first unvisited neighbor
      int neighbor = neighbors[0];

      // Add edge between current node and neighbor in the tree
      tree_adj_matrix(current_node, neighbor) = 1;
      tree_adj_matrix(neighbor, current_node) = 1;  // Ensure symmetry

      // Mark the neighbor as visited and push it to the stack
      visited[neighbor] = true;
      stack.push_back(current_node);  // Revisit the current node for other neighbors
      stack.push_back(neighbor);      // Visit the new neighbor next

      edges_added++;
    }
  }

  return tree_adj_matrix;  // Return the adjacency matrix of the random spanning tree
}

struct Edge {
  double weight;
  int u, v;
  bool operator>(const Edge& other) const { return weight > other.weight; }
};

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
Rcpp::List generate_random_spanning_tree_and_probability(const arma::mat& weight_matrix) {
  int n = weight_matrix.n_rows;
  arma::mat tree_adj_matrix = arma::zeros<arma::mat>(n, n);  // Initialize spanning tree adjacency matrix
  double tree_weight_product = 1.0;  // Product of the edge weights in the spanning tree

  std::vector<bool> in_tree(n, false);  // Track if a node is in the tree
  std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> edge_queue;  // Priority queue for edges

  // Random number generator
  std::random_device rd;
  std::mt19937 gen(rd());

  // Start from a random node
  std::uniform_int_distribution<> distrib(0, n-1);
  int start_node = distrib(gen);
  in_tree[start_node] = true;

  // Push all edges connected to the start node into the queue
  for (int j = 0; j < n; ++j) {
    if (weight_matrix(start_node, j) > 0) {  // Only consider positive-weight edges
      edge_queue.push({weight_matrix(start_node, j), start_node, j});
    }
  }

  // Keep track of how many edges we've added to the tree
  int edges_added = 0;

  // Continue until the tree has n-1 edges
  while (!edge_queue.empty() && edges_added < n - 1) {
    // Extract the smallest edge from the queue
    Edge current_edge = edge_queue.top();
    edge_queue.pop();

    // Check if the edge forms a cycle (i.e., both nodes are already in the tree)
    if (in_tree[current_edge.u] && in_tree[current_edge.v]) {
      continue;  // Skip this edge if it forms a cycle
    }

    // Add the edge to the spanning tree
    tree_adj_matrix(current_edge.u, current_edge.v) = 1;
    tree_adj_matrix(current_edge.v, current_edge.u) = 1;
    edges_added++;

    // Multiply the weight of the current edge to the total weight product
    tree_weight_product *= current_edge.weight;

    // Add the new node to the tree
    int new_node = in_tree[current_edge.u] ? current_edge.v : current_edge.u;
    in_tree[new_node] = true;

    // Push all edges from the new node into the queue
    for (int j = 0; j < n; ++j) {
      if (!in_tree[j] && weight_matrix(new_node, j) > 0) {
        edge_queue.push({weight_matrix(new_node, j), new_node, j});
      }
    }
  }

  // The probability of the tree is proportional to the product of the edge weights
  double probability = tree_weight_product;

  return Rcpp::List::create(
    Rcpp::Named("adjacency_matrix") = tree_adj_matrix,
    Rcpp::Named("tree_weight_product") = tree_weight_product,
    Rcpp::Named("probability") = probability
  );
}

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
Rcpp::List generate_random_spanning_tree_and_log_probability(const arma::mat& weight_matrix) {
  int n = weight_matrix.n_rows;
  arma::mat tree_adj_matrix = arma::zeros<arma::mat>(n, n);  // Initialize spanning tree adjacency matrix
  double tree_weight_product = 1.0;  // Product of the edge weights in the spanning tree

  std::vector<bool> in_tree(n, false);  // Track if a node is in the tree
  std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> edge_queue;  // Priority queue for edges

  // Random number generator
  std::random_device rd;
  std::mt19937 gen(rd());

  // Start from a random node
  std::uniform_int_distribution<> distrib(0, n-1);
  int start_node = distrib(gen);
  in_tree[start_node] = true;

  // Push all edges connected to the start node into the queue
  for (int j = 0; j < n; ++j) {
    if (weight_matrix(start_node, j) > 0) {  // Only consider positive-weight edges
      edge_queue.push({weight_matrix(start_node, j), start_node, j});
    }
  }

  // Keep track of how many edges we've added to the tree
  int edges_added = 0;

  // Continue until the tree has n-1 edges
  while (!edge_queue.empty() && edges_added < n - 1) {
    // Extract the smallest edge from the queue
    Edge current_edge = edge_queue.top();
    edge_queue.pop();

    // Check if the edge forms a cycle (i.e., both nodes are already in the tree)
    if (in_tree[current_edge.u] && in_tree[current_edge.v]) {
      continue;  // Skip this edge if it forms a cycle
    }

    // Add the edge to the spanning tree
    tree_adj_matrix(current_edge.u, current_edge.v) = 1;
    tree_adj_matrix(current_edge.v, current_edge.u) = 1;
    edges_added++;

    // Multiply the weight of the current edge to the total weight product
    tree_weight_product += log(current_edge.weight);

    // Add the new node to the tree
    int new_node = in_tree[current_edge.u] ? current_edge.v : current_edge.u;
    in_tree[new_node] = true;

    // Push all edges from the new node into the queue
    for (int j = 0; j < n; ++j) {
      if (!in_tree[j] && weight_matrix(new_node, j) > 0) {
        edge_queue.push({weight_matrix(new_node, j), new_node, j});
      }
    }
  }

  // The probability of the tree is proportional to the product of the edge weights
  double probability = tree_weight_product;

  return Rcpp::List::create(
    Rcpp::Named("adjacency_matrix") = tree_adj_matrix,
    Rcpp::Named("tree_weight_product") = tree_weight_product
    //Rcpp::Named("probability") = probability
  );
}

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double jaccard_similarity(const arma::mat& adj1, const arma::mat& adj2) {
  if (adj1.n_rows != adj2.n_rows || adj1.n_cols != adj2.n_cols) {
    stop("Adjacency matrices must have the same dimensions.");
  }

  // Flatten the adjacency matrices into vectors (excluding diagonal elements for undirected graphs)
  arma::uvec edges1 = arma::find(adj1);  // Get indices of non-zero elements (edges) in adj1
  arma::uvec edges2 = arma::find(adj2);  // Get indices of non-zero elements (edges) in adj2

  // Calculate intersection (common edges)
  arma::uvec intersection = arma::intersect(edges1, edges2);

  // Calculate union (total distinct edges)
  arma::uvec union_edges = arma::unique(arma::join_cols(edges1, edges2));

  // Jaccard similarity is the size of the intersection divided by the size of the union
  double jaccard_sim = (double)intersection.n_elem / (double)union_edges.n_elem;

  return jaccard_sim;
}

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List bayestrandgrap_cpp(const arma::vec& y, const arma::mat& xobs, const arma::vec& thetaini,
                        int iter, int burn, int thin, const arma::mat& weight_mat,
                        const arma::vec& aprior, const arma::vec& bprior,
                        double divproprho = 18, double divproppsi = 20) {

  int n = y.n_elem;
  int p = xobs.n_cols;

  // Initialize vectors and lists
  arma::vec rhoF(iter, fill::zeros);
  arma::vec psiF(iter, fill::zeros);
  std::vector<arma::mat> WF(iter);
  std::vector<double> logpw(iter);

  rhoF(0) = thetaini(0);
  psiF(0) = thetaini(1);

  // Initial spanning tree
  List init_tree = generate_random_spanning_tree_and_log_probability(weight_mat);
  arma::mat adjsamini = as<arma::mat>(init_tree["adjacency_matrix"]);
  double init_logpw = as<double>(init_tree["tree_weight_product"]);

  WF[0] = adjsamini;
  logpw[0] = init_logpw;

  int count = 0;

  // MCMC loop
  for (int i = 1; i < iter; ++i) {
    List cand_tree = generate_random_spanning_tree_and_log_probability(weight_mat);
    arma::mat adjsamcand = as<arma::mat>(cand_tree["adjacency_matrix"]);
    double logpwcand = as<double>(cand_tree["tree_weight_product"]);

    arma::mat adjsamlast = WF[i-1];
    double logpwlast = logpw[i-1];

    // Proposal for rho and psi
    arma::vec thetalast = {rhoF[i-1], psiF[i-1]};
    arma::vec muprop = {std::min(0.99, thetalast[0]), std::min(0.99, thetalast[1])};
    arma::vec sigmaprop = {muprop[0] * (1 - muprop[0]) / divproprho,
                           muprop[1] * (1 - muprop[1]) / divproppsi};

    double omega = arma::randu() * (0.99 - 1.0) + 1.0;
    double theta1cand = std::min(0.99, rbeta_rep_cpp(omega * muprop[0], sigmaprop[0]));
    double theta2cand = std::min(0.99, rbeta_rep_cpp(omega * muprop[1], sigmaprop[1]));
    arma::vec thetacand = {theta1cand, theta2cand};

    // Posterior ratio (log-scale)
    double lognum_aux = log(posteriorvar3_cpp(thetacand, adjsamcand, xobs, y, aprior, bprior)) + logpwcand;
    double logden_aux = log(posteriorvar3_cpp(thetalast, adjsamlast, xobs, y, aprior, bprior)) + logpwlast;

    double lognum = lognum_aux + log(dbeta_rep_cpp(thetalast[0],omega*muprop[0],sigmaprop[0])) + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + logpwlast;
    double logden = logden_aux + log(dbeta_rep_cpp(thetacand[0],omega*muprop[0],sigmaprop[0])) + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + logpwcand;


    // Metropolis-Hastings step
    if (std::log(arma::randu()) < lognum - logden) {
      rhoF[i] = theta1cand;
      psiF[i] = theta2cand;
      WF[i] = adjsamcand;
      logpw[i] = logpwcand;
      ++count;
    } else {
      rhoF[i] = rhoF[i-1];
      psiF[i] = psiF[i-1];
      WF[i] = adjsamlast;
      logpw[i] = logpwlast;
    }
  }

  // Post-processing
  arma::vec rhoburn = rhoF.subvec(burn, iter-1);
  arma::vec rhoval = rhoburn.subvec(0, rhoburn.n_elem / thin - 1);

  arma::vec psiburn = psiF.subvec(burn, iter-1);
  arma::vec psival = psiburn.subvec(0, psiburn.n_elem / thin - 1);

  // You would complete the rest of the logic similar to the original R code here
  // Sampling betaFval, sigmaFval etc.

  return List::create(Named("rhoF") = rhoval,
                      Named("psiF") = psival,
                      Named("probacc") = count);
}


// Precompute neighbors and probabilities for each node
std::vector<std::vector<int>> precompute_neighbors(const arma::mat& weight_matrix, std::vector<std::vector<double>>& probabilities) {
  int n = weight_matrix.n_rows;
  std::vector<std::vector<int>> neighbors(n);

  for (int i = 0; i < n; ++i) {
    double total_weight = 0.0;
    for (int j = 0; j < n; ++j) {
      if (weight_matrix(i, j) > 0) {
        neighbors[i].push_back(j);
        probabilities[i].push_back(weight_matrix(i, j));
        total_weight += weight_matrix(i, j); // Sum the weights of the row
      }
    }
    // Normalize probabilities for node i
    for (auto& p : probabilities[i]) {
      p /= total_weight;  // Normalize to make them probabilities
    }
  }

  return neighbors;
}

// A helper function to choose the next node in the biased random walk (Aldous-Broder)
int choose_next_node_aldous_broder(int current_node, const std::vector<std::vector<int>>& neighbors,
                                   const std::vector<std::vector<double>>& probabilities, std::mt19937& gen) {
  // Set up a discrete distribution for choosing the next node
  std::discrete_distribution<int> dist(probabilities[current_node].begin(), probabilities[current_node].end());
  return neighbors[current_node][dist(gen)];
}

//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
Rcpp::List generate_weighted_spanning_tree_aldous_broder_optimized(const arma::mat& weight_matrix) {
  int n = weight_matrix.n_rows;

  // Initialize the spanning tree adjacency matrix
  arma::sp_mat tree_adj_matrix(n, n);  // Use a sparse matrix to store the tree

  // Track which nodes are already part of the spanning tree
  std::vector<bool> in_tree(n, false);

  // Precompute neighbors and probabilities
  std::vector<std::vector<double>> probabilities(n);
  std::vector<std::vector<int>> neighbors = precompute_neighbors(weight_matrix, probabilities);

  // Random number generator
  std::random_device rd;
  std::mt19937 gen(rd());

  // Start the random walk from a random node
  std::uniform_int_distribution<> distrib(0, n - 1);
  int current_node = distrib(gen);
  in_tree[current_node] = true;

  // Probability of the spanning tree (logarithmic)
  double log_tree_weight = 0.0;

  // Total number of edges added to the tree
  int edges_added = 0;

  while (edges_added < n - 1) {
    // Choose the next node in the random walk based on precomputed neighbors and probabilities
    int next_node = choose_next_node_aldous_broder(current_node, neighbors, probabilities, gen);

    // If the next node is not yet in the tree, add the edge and update the tree
    if (!in_tree[next_node]) {
      tree_adj_matrix(current_node, next_node) = 1;
      tree_adj_matrix(next_node, current_node) = 1;
      in_tree[next_node] = true;

      // Compute the log probability using the ratio of the edge weight to the sum of the weights in the current node's row
      double row_sum = 0.0;
      for (int j = 0; j < n; ++j) {
        row_sum += weight_matrix(current_node, j); // Sum of weights in the current node's row
      }
      double edge_weight = weight_matrix(current_node, next_node);
      log_tree_weight += std::log(edge_weight / row_sum);  // Logarithm of the ratio

      edges_added++;
    }

    // Move to the next node in the random walk
    current_node = next_node;
  }

  // The final probability of the tree is exp(log_tree_weight)
  return Rcpp::List::create(
    Rcpp::Named("adjacency_matrix") = tree_adj_matrix,
    Rcpp::Named("log_probability") = log_tree_weight,
    Rcpp::Named("probability") = std::exp(log_tree_weight)
  );
}


// Precompute neighbors and probabilities for each node, including row sums
std::vector<std::vector<int>> precompute_neighbors_2(const arma::mat& weight_matrix,
                                                     std::vector<std::vector<double>>& probabilities,
                                                     std::vector<double>& row_sums) {
  int n = weight_matrix.n_rows;
  std::vector<std::vector<int>> neighbors(n);

  for (int i = 0; i < n; ++i) {
    double total_weight = 0.0;
    for (int j = 0; j < n; ++j) {
      if (weight_matrix(i, j) > 0) {
        neighbors[i].push_back(j);
        double probability = weight_matrix(i, j);
        probabilities[i].push_back(probability);
        total_weight += probability;
      }
    }
    row_sums[i] = total_weight;

    // Normalize probabilities in-place for node i
    for (double &p : probabilities[i]) {
      p /= total_weight;
    }
  }
  return neighbors;
}

// A helper function to choose the next node in the biased random walk (Aldous-Broder)
int choose_next_node_aldous_broder_2(int current_node, const std::vector<std::vector<int>>& neighbors,
                                     const std::vector<std::vector<double>>& probabilities, std::mt19937& gen) {
  // Set up a discrete distribution for choosing the next node
  std::discrete_distribution<int> dist(probabilities[current_node].begin(), probabilities[current_node].end());
  return neighbors[current_node][dist(gen)];
}

//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
 // [[Rcpp::export]]
 Rcpp::List generate_weighted_spanning_tree_aldous_broder_optimized_2(const arma::mat& weight_matrix) {
   int n = weight_matrix.n_rows;

   // Initialize the spanning tree adjacency matrix as sparse
   arma::sp_mat tree_adj_matrix(n, n);

   // Track which nodes are already part of the spanning tree
   std::vector<bool> in_tree(n, false);

   // Precompute neighbors, probabilities, and row sums
   std::vector<std::vector<double>> probabilities(n);
   std::vector<double> row_sums(n);
   std::vector<std::vector<int>> neighbors = precompute_neighbors_2(weight_matrix, probabilities, row_sums);

   // Random number generator
   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_int_distribution<> distrib(0, n - 1);

   // Start the random walk from a random node
   int current_node = distrib(gen);
   in_tree[current_node] = true;

   // Probability of the spanning tree (logarithmic)
   double log_tree_weight = 0.0;

   // Total number of edges added to the tree
   int edges_added = 0;

   while (edges_added < n - 1) {
     // Choose the next node in the random walk
     int next_node = choose_next_node_aldous_broder_2(current_node, neighbors, probabilities, gen);

     // If the next node is not yet in the tree, add the edge and update the tree
     if (!in_tree[next_node]) {
       tree_adj_matrix(current_node, next_node) = 1;
       tree_adj_matrix(next_node, current_node) = 1;
       in_tree[next_node] = true;

       // Compute the log probability using the precomputed row sum
       double edge_weight = weight_matrix(current_node, next_node);
       log_tree_weight += std::log(edge_weight / row_sums[current_node]);

       edges_added++;
     }

     // Move to the next node in the random walk
     current_node = next_node;
   }

   // The final probability of the tree is exp(log_tree_weight)
   return Rcpp::List::create(
     Rcpp::Named("adjacency_matrix") = tree_adj_matrix,
     Rcpp::Named("log_probability") = log_tree_weight//,
    // Rcpp::Named("probability") = std::exp(log_tree_weight)
   );
 }

//' @useDynLib auxiliarcpp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::vec loglikespatemp_cpp(List listres, const arma::mat& Xmat, const arma::vec& y, int lags) {
  int p = Xmat.n_cols;
  int N = y.n_elem;

  arma::vec theta = listres["thetaFval"];
  arma::mat adjmatinf = listres["adjmatinf"];

  arma::vec beta = theta.subvec(0, p - 1);
  double sigma2 = theta(p);
  double rho = theta(p + 1);
  double psi = theta(p + 2);
  double phi = theta(p + 3);
  arma::vec mufix = Xmat * beta;

  // Calculate covariance matrix
  arma::mat varcova = sigma2 * spatimecovar_2_cpp(lags, adjmatinf, rho, psi, phi); // Replace 0 with phi if needed

  // Inverse covariance matrix
  arma::mat invvarcova = arma::inv(varcova);

  // Calculate sigmamu
  arma::vec sigmamu = invvarcova * (y - mufix);

  arma::vec loglik(N);
  for (int i = 0; i < N; ++i) {
    double sigmai = 1.0 / std::sqrt(invvarcova(i, i));
    double mui = y(i) - (sigmamu(i) * sigmai * sigmai);
    loglik(i) = R::dnorm(y(i), mui, sigmai, true);
  }

  return loglik;
}



// ---- utilities: force clean R numeric vectors/matrices ----
static inline NumericVector to_nv(const arma::vec& v) {
  NumericVector out(v.n_elem);
  std::copy(v.begin(), v.end(), out.begin());
  out.attr("names") = R_NilValue;
  return out;
}
static inline NumericMatrix to_nm(const arma::mat& M) {
  NumericMatrix out(M.n_rows, M.n_cols);
  std::copy(M.begin(), M.end(), out.begin()); // column-major OK
  out.attr("dimnames") = R_NilValue;
  return out;
}

// ---- rectangular MVN log-prob via mvtnorm (robust coercions) ----
static double log_rect_mvn_prob(const arma::vec& lower,
                                const arma::vec& upper,
                                const arma::vec& mean,
                                const arma::mat& Sigma) {
  static Environment mvtnorm = Environment::namespace_env("mvtnorm");
  static Function pmvnorm    = mvtnorm["pmvnorm"];
  static Function GenzBretz  = mvtnorm["GenzBretz"];

  NumericVector lower_nv = to_nv(lower);
  NumericVector upper_nv = to_nv(upper);
  NumericVector mean_nv  = to_nv(mean);
  NumericMatrix Sigma_nm = to_nm(Sigma);

  List alg = GenzBretz(_["maxpts"] = 10000, _["abseps"] = 1e-6, _["releps"] = 0.0);

  NumericVector prob_sexp = pmvnorm(_["lower"]= lower_nv,
                                    _["upper"]= upper_nv,
                                    _["mean"] = mean_nv,
                                    _["sigma"]= Sigma_nm,
                                    _["algorithm"]= alg);
  if (prob_sexp.size() != 1) stop("pmvnorm did not return a scalar.");
  double prob = prob_sexp[0];
  if (!std::isfinite(prob) || prob <= 0.0)
    return -std::numeric_limits<double>::infinity();
  return std::log(prob);
}

// ---- main: build covariance by calling a function in a given package namespace ----
// pass pkg="auxiliarcpp", fun_name="spatimecovar2" (or the exact symbol).
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
// [[Rcpp::export]]
double log_likelihood_ar1_fast_pkg_cpp(const arma::mat& x,           // n x p
                                       const arma::vec& y,           // n
                                       const arma::ivec& cc,         // 0 obs, 1 cens
                                       int lag,
                                       const arma::mat& adjmatinf,
                                       const arma::vec& theta,       // [beta(1:p), sigma2, rho, psi, phi]
                                       const arma::vec& lower,       // #cens
                                       const arma::vec& upper,       // #cens
                                       std::string pkg,
                                       std::string fun_name = "spatimecovar_2_cpp") {

  const int p = x.n_cols;
  if ((int)theta.n_elem < p + 4)
    stop("theta must contain beta(1:p), sigma2, rho, psi, phi");

  arma::vec  beta   = theta.subvec(0, p - 1);
  double     sigma2 = theta[p];
  double     rho    = theta[p + 1];
  double     psi    = theta[p + 2];
  double     phi    = theta[p + 3];

  arma::vec eta = x * beta;

  arma::uvec idx_obs  = arma::find(cc == 0);
  arma::uvec idx_cens = arma::find(cc == 1);
  const std::size_t n_obs  = idx_obs.n_elem;
  const std::size_t n_cens = idx_cens.n_elem;

  if ((std::size_t)lower.n_elem != n_cens || (std::size_t)upper.n_elem != n_cens)
    stop("lower/upper must have length equal to the number of censored entries.");

  // get covariance builder from package namespace
  Environment penv = Environment::namespace_env(pkg);
  if (!penv.exists(fun_name)) {
    stop("Function '%s' not found in package namespace '%s'.",
         fun_name.c_str(), pkg.c_str());
  }
  Function covfun = penv[fun_name];

  // Build covariance and scale by sigma^2
  arma::mat cov_full = as<arma::mat>(covfun(_["lag"]= lag,
                                            _["adjmatinf"]= wrap(adjmatinf),
                                            _["rho"]= rho, _["psi"]= psi, _["phi"]= phi));
  cov_full *= sigma2;

  // no observed? (all censored or none)
  if (n_obs == 0) {
    if (n_cens == 0) return 0.0;
    arma::vec mu_c    = eta.elem(idx_cens);
    arma::mat Sigma_c = arma::symmatu(cov_full.submat(idx_cens, idx_cens));
    return log_rect_mvn_prob(lower, upper, mu_c, Sigma_c);
  }

  // observed block via single Cholesky
  arma::vec y_obs   = y.elem(idx_obs);
  arma::vec eta_obs = eta.elem(idx_obs);
  arma::mat Sigma_oo = arma::symmatu(cov_full.submat(idx_obs, idx_obs));

  arma::mat L;
  if (!arma::chol(L, Sigma_oo, "lower")) {
    double jitter = 1e-10 * arma::trace(Sigma_oo) /
      std::max<std::size_t>((std::size_t)1, Sigma_oo.n_rows);
    if (!arma::chol(L, Sigma_oo + jitter * arma::eye(Sigma_oo.n_rows, Sigma_oo.n_cols), "lower")) {
      stop("Cholesky failed for Sigma_oo (even after jitter).");
    }
  }

  arma::vec r = y_obs - eta_obs;
  arma::vec z = arma::solve(arma::trimatl(L), r);

  double logdet = 2.0 * arma::sum(arma::log(L.diag()));
  double loglik = -0.5 * ( static_cast<double>(n_obs) * std::log(2.0 * arma::datum::pi)
                             + logdet + arma::dot(z, z) );

  // censored conditional
  if (n_cens > 0) {
    arma::vec eta_cen = eta.elem(idx_cens);
    arma::mat Sigma_oc = cov_full.submat(idx_obs,  idx_cens);
    arma::mat Sigma_cc = arma::symmatu(cov_full.submat(idx_cens, idx_cens));

    arma::mat K = arma::solve(arma::trimatl(L), Sigma_oc);   // L * K = Sigma_oc
    arma::vec mu_c = eta_cen + K.t() * z;
    arma::mat Sigma_c = arma::symmatu(Sigma_cc - K.t() * K);

    loglik += log_rect_mvn_prob(lower, upper, mu_c, Sigma_c);
  }

  return loglik;
}



//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
// [[Rcpp::export]]
double log_likelihood_ar2_fast_pkg_cpp(const arma::mat& x,           // n x p
                                       const arma::vec& y,           // n
                                       const arma::ivec& cc,         // 0 obs, 1 cens
                                       int lag,
                                       const arma::mat& adjmatinf,
                                       const arma::vec& theta,       // [beta(1:p), sigma2, rho, psi, phi1, phi2]
                                       const arma::vec& lower,       // #cens
                                       const arma::vec& upper,       // #cens
                                       std::string pkg,
                                       std::string fun_name = "spatimecovar_ar2_cpp") {

  const int p = x.n_cols;
  if ((int)theta.n_elem < p + 5)
    stop("theta must contain beta(1:p), sigma2, rho, psi, phi1, phi2");

  arma::vec  beta   = theta.subvec(0, p - 1);
  double     sigma2 = theta[p];
  double     rho    = theta[p + 1];
  double     psi    = theta[p + 2];
  double     phi1   = theta[p + 3];
  double     phi2   = theta[p + 4];

  arma::vec eta = x * beta;

  arma::uvec idx_obs  = arma::find(cc == 0);
  arma::uvec idx_cens = arma::find(cc == 1);
  const std::size_t n_obs  = idx_obs.n_elem;
  const std::size_t n_cens = idx_cens.n_elem;

  if ((std::size_t)lower.n_elem != n_cens || (std::size_t)upper.n_elem != n_cens)
    stop("lower/upper must have length equal to the number of censored entries.");

  // get covariance builder (spatimecovar_ar2_cpp) from package namespace
  Environment penv = Environment::namespace_env(pkg);
  if (!penv.exists(fun_name)) {
    stop("Function '%s' not found in package namespace '%s'.",
         fun_name.c_str(), pkg.c_str());
  }
  Function covfun = penv[fun_name];

  // Build covariance (returns timespa + psi*I) and scale by sigma^2
  arma::mat cov_full = as<arma::mat>(covfun(_["lag"]= lag,
                                            _["adjmatinf"]= wrap(adjmatinf),
                                            _["rho"]= rho, _["psi"]= psi,
                                            _["phi1"]= phi1, _["phi2"]= phi2));
  cov_full *= sigma2;

  // all censored or none
  if (n_obs == 0) {
    if (n_cens == 0) return 0.0;
    arma::vec mu_c    = eta.elem(idx_cens);
    arma::mat Sigma_c = arma::symmatu(cov_full.submat(idx_cens, idx_cens));
    return log_rect_mvn_prob(lower, upper, mu_c, Sigma_c);
  }

  // observed block via single Cholesky
  arma::vec y_obs    = y.elem(idx_obs);
  arma::vec eta_obs  = eta.elem(idx_obs);
  arma::mat Sigma_oo = arma::symmatu(cov_full.submat(idx_obs, idx_obs));

  arma::mat L;
  if (!arma::chol(L, Sigma_oo, "lower")) {
    double jitter = 1e-10 * arma::trace(Sigma_oo) /
      std::max<std::size_t>((std::size_t)1, Sigma_oo.n_rows);
    if (!arma::chol(L, Sigma_oo + jitter * arma::eye(Sigma_oo.n_rows, Sigma_oo.n_cols), "lower")) {
      stop("Cholesky failed for Sigma_oo (even after jitter).");
    }
  }

  arma::vec r = y_obs - eta_obs;
  arma::vec z = arma::solve(arma::trimatl(L), r);

  double logdet = 2.0 * arma::sum(arma::log(L.diag()));
  double loglik = -0.5 * ( static_cast<double>(n_obs) * std::log(2.0 * arma::datum::pi)
                             + logdet + arma::dot(z, z) );

  // censored conditional contribution
  if (n_cens > 0) {
    arma::vec eta_c   = eta.elem(idx_cens);
    arma::mat Sigma_oc = cov_full.submat(idx_obs,  idx_cens);
    arma::mat Sigma_cc = arma::symmatu(cov_full.submat(idx_cens, idx_cens));

    arma::mat K = arma::solve(arma::trimatl(L), Sigma_oc);   // L * K = Sigma_oc
    arma::vec mu_c = eta_c + K.t() * z;
    arma::mat Sigma_c = arma::symmatu(Sigma_cc - K.t() * K);

    loglik += log_rect_mvn_prob(lower, upper, mu_c, Sigma_c);
  }

  return loglik;
}


//// CAR MODEL

/// taken FUNCTION TGMRF package

arma::mat buildQST_cpp(const arma::mat& Ws,
                       const arma::mat& Wt,
                       const double     rho_s,
                       const double     rho_t,
                       const double     rho_st) {

  const int nReg = Ws.n_rows;
  const int nVar = Wt.n_rows;

  arma::mat Inreg = arma::eye<arma::mat>(nReg, nReg);
  arma::mat Invar = arma::eye<arma::mat>(nVar, nVar);

  arma::mat kronAux1 = arma::kron(Ws,   Invar);
  arma::mat kronAux2 = arma::kron(Inreg,Wt);
  arma::mat kronAux3 = arma::kron(Ws,   Wt);

  arma::mat Q = -(rho_s * kronAux1 + rho_t * kronAux2 + rho_st * kronAux3);

  // Column-wise cumulative sums; last row of each column equals its column sum
  arma::mat D = arma::cumsum(kronAux1 + kronAux2 + kronAux3, 0);

  const int nn = nReg * nVar;
  for (int i = 0; i < nn; ++i)
    Q(i, i) = D(nn - 1, i);  // set diagonal to degree (column sums)

  return Q;
}

/*** Space–time covariance-like matrix (Q + psi I).
 Keep ONE exported signature that you actually use. ***/
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
// [[Rcpp::export]]
arma::mat spatimecovarcar_cpp(const arma::mat&   adjmatinf,
                              const arma::mat&   Wt,
                              const double       rho_s,
                              const double       rho_t,
                              const double       rho_st,
                              const double       psi) {
  //(void)lag; // remove if you later use it

  arma::mat Q = buildQST_cpp(adjmatinf, Wt, rho_s, rho_t, rho_st);
  arma::mat I = arma::eye<arma::mat>(Q.n_rows, Q.n_cols);
  return arma::inv(Q) + psi * I;
}


//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::mat spatimecovarcar_zero_cpp(const arma::mat&   adjmatinf,
                               const arma::mat&   Wt,
                               const double       rho_s,
                               const double       rho_t,
                               const double       rho_st) {
   //(void)lag; // remove if you later use it

   arma::mat Q = buildQST_cpp(adjmatinf, Wt, rho_s, rho_t, rho_st);
   arma::mat I = arma::eye<arma::mat>(Q.n_rows, Q.n_cols);
   return arma::inv(Q);
 }








Rcpp::List posterior_spatempcenscar2_cpp(const arma::vec& theta,
                                         const arma::mat& adjmatinf,
                                         const arma::mat&   Wt,
                                         const arma::mat& X,
                                         const arma::vec& y,
                                         const arma::vec& a,
                                         const arma::vec& b,
                                         const arma::vec& alim,
                                         const arma::vec& blim);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 Rcpp::List posterior_spatempcenscar2_cpp(const arma::vec& theta,
                                          const arma::mat& adjmatinf,
                                          const arma::mat&   Wt,
                                          const arma::mat& X,
                                          const arma::vec& y,
                                          const arma::vec& a,
                                          const arma::vec& b,
                                          const arma::vec& alim,
                                          const arma::vec& blim) {
   double rho_s = theta[0];
   double rho_t = theta[1];
   double rho_st = theta[2];
   double psi = theta[3];

   double arho_s = a[0];
   double arho_t = a[1];
   double arho_st = a[2];
   double apsi = a[3];

   double brho_s = b[0];
   double brho_t = b[1];
   double brho_st = b[2];
   double bpsi = b[3];

   double alimrho_s = alim[0];
   double alimrho_t = alim[1];
   double alimrho_st = alim[2];
   double alimpsi = alim[3];

   double blimrho_s = blim[0];
   double blimrho_t = blim[1];
   double blimrho_st = blim[2];
   double blimpsi = blim[3];


   // Compute prior densities for rho and psi
   double priorrho_s =  dbeta_ab(rho_s, arho_s, brho_s, alimrho_s, blimrho_s,false);
   double priorrho_t =  dbeta_ab(rho_t, arho_t, brho_t, alimrho_t, blimrho_t,false);
   double priorrho_st = dbeta_ab(rho_st, arho_st, brho_st, alimrho_st, blimrho_st,false);
   double priorpsi =    dbeta_ab(psi, apsi, bpsi, alimpsi, blimpsi, false);


   // Compute variance-covariance matrix
   arma::mat varcovdagar = spatimecovarcar_cpp(adjmatinf,Wt, rho_s, rho_t, rho_st, psi);
   arma::mat invvarcovdagar = arma::inv(varcovdagar);

   arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
   arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


   int n = X.n_rows;
   int p = X.n_cols;

   // Compute residuals and variance
   arma::vec residuals = y - X * betaest;
   double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

   // Determinants for posterior calculation
   double detvarcovinv = arma::det(invvarcovdagar);
   if (detvarcovinv == R_PosInf) {
     detvarcovinv = 1e308;
   }
   if (detvarcovinv == 0) {
     detvarcovinv = 1e-300;
   }
   double log_detvarcovinv = log(detvarcovinv);
   double detvarbeta = arma::det(varbeta);
   double log_detvarbeta = log(detvarbeta);
   double s2np = std::pow(s2, -(n - p));
   double log_s2np = -(n-p)*log(s2);

   // Return the posterior variance
   double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho_s*priorrho_t*priorrho_st* priorpsi;
   double posterior_log = 0.5*(log_detvarcovinv + log_detvarbeta +log_s2np) + log(priorrho_s) +log(priorrho_t) + log(priorrho_st)  +log(priorpsi);
   return Rcpp::List::create(
     Rcpp::Named("priorrho_s") = priorrho_s,
     Rcpp::Named("priorrho_t") = priorrho_t,
     Rcpp::Named("priorrho_st") = priorrho_st,
     Rcpp::Named("priorpsi") = priorpsi,
     Rcpp::Named("varcovspatemp") = varcovdagar,
     Rcpp::Named("varbeta") = varbeta,
     Rcpp::Named("betaest") = betaest,
     Rcpp::Named("S2") = s2,
     Rcpp::Named("posterior") = posterior,
     Rcpp::Named("posterior_log") = posterior_log
   );
   //List result;
   //result["detvcovinv"] = detvarcovinv;
   //result["detvarbeta"] = detvarbeta;
   //result["s2expo"] = std::pow(s2, -(n - p));
   //result["posterior"] = posterior;
   //result["priorrho"] = priorrho;
   //result["priorpsi"] = priorpsi;
   //result["priorphi"] = priorphi;
   //return result;
 }



Rcpp::List posterior_spatempcenscar_zero_cpp(const arma::vec& theta,
                                         const arma::mat& adjmatinf,
                                         const arma::mat&   Wt,
                                         const arma::mat& X,
                                         const arma::vec& y,
                                         const arma::vec& a,
                                         const arma::vec& b,
                                         const arma::vec& alim,
                                         const arma::vec& blim);
//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 Rcpp::List posterior_spatempcenscar_zero_cpp(const arma::vec& theta,
                                          const arma::mat& adjmatinf,
                                          const arma::mat&   Wt,
                                          const arma::mat& X,
                                          const arma::vec& y,
                                          const arma::vec& a,
                                          const arma::vec& b,
                                          const arma::vec& alim,
                                          const arma::vec& blim) {
   double rho_s = theta[0];
   double rho_t = theta[1];
   double rho_st = theta[2];

   double arho_s = a[0];
   double arho_t = a[1];
   double arho_st = a[2];

   double brho_s = b[0];
   double brho_t = b[1];
   double brho_st = b[2];

   double alimrho_s = alim[0];
   double alimrho_t = alim[1];
   double alimrho_st = alim[2];

   double blimrho_s = blim[0];
   double blimrho_t = blim[1];
   double blimrho_st = blim[2];


   // Compute prior densities for rho and psi
   double priorrho_s =  dbeta_ab(rho_s, arho_s, brho_s, alimrho_s, blimrho_s,false);
   double priorrho_t =  dbeta_ab(rho_t, arho_t, brho_t, alimrho_t, blimrho_t,false);
   double priorrho_st = dbeta_ab(rho_st, arho_st, brho_st, alimrho_st, blimrho_st,false);


   // Compute variance-covariance matrix
   arma::mat varcovdagar = spatimecovarcar_zero_cpp(adjmatinf,Wt, rho_s, rho_t, rho_st);
   arma::mat invvarcovdagar = arma::inv(varcovdagar);

   arma::mat varbeta = arma::inv(X.t() * invvarcovdagar * X);
   arma::vec betaest = varbeta * X.t() * invvarcovdagar * y;


   int n = X.n_rows;
   int p = X.n_cols;

   // Compute residuals and variance
   arma::vec residuals = y - X * betaest;
   double s2 = arma::as_scalar(residuals.t() * invvarcovdagar * residuals) / (n - p);

   // Determinants for posterior calculation
   double detvarcovinv = arma::det(invvarcovdagar);
   if (detvarcovinv == R_PosInf) {
     detvarcovinv = 1e308;
   }
   if (detvarcovinv == 0) {
     detvarcovinv = 1e-300;
   }
   double log_detvarcovinv = log(detvarcovinv);
   double detvarbeta = arma::det(varbeta);
   double log_detvarbeta = log(detvarbeta);
   double s2np = std::pow(s2, -(n - p));
   double log_s2np = -(n-p)*log(s2);

   // Return the posterior variance
   double posterior = std::sqrt(detvarcovinv*detvarbeta * std::pow(s2, -(n - p))) * priorrho_s*priorrho_t*priorrho_st;
   double posterior_log = 0.5*(log_detvarcovinv + log_detvarbeta +log_s2np) + log(priorrho_s) +log(priorrho_t) + log(priorrho_st);
   return Rcpp::List::create(
     Rcpp::Named("priorrho_s") = priorrho_s,
     Rcpp::Named("priorrho_t") = priorrho_t,
     Rcpp::Named("priorrho_st") = priorrho_st,
     Rcpp::Named("varcovspatemp") = varcovdagar,
     Rcpp::Named("varbeta") = varbeta,
     Rcpp::Named("betaest") = betaest,
     Rcpp::Named("S2") = s2,
     Rcpp::Named("posterior") = posterior,
     Rcpp::Named("posterior_log") = posterior_log
   );
   //List result;
   //result["detvcovinv"] = detvarcovinv;
   //result["detvarbeta"] = detvarbeta;
   //result["s2expo"] = std::pow(s2, -(n - p));
   //result["posterior"] = posterior;
   //result["priorrho"] = priorrho;
   //result["priorpsi"] = priorpsi;
   //result["priorphi"] = priorphi;
   //return result;
 }


// simple AR(1)-like time adjacency: 1 on |i-j|=1
static inline arma::mat make_Wt(const int lag) {
  arma::mat Wt(lag, lag, fill::zeros);
  for (int i = 0; i + 1 < lag; ++i) {
    Wt(i, i + 1) = 1.0;
    Wt(i + 1, i) = 1.0;
  }
  return Wt;
}

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
// [[Rcpp::export]]
double log_likelihood_car_fast_cpp(const arma::mat& x,
                                   const arma::vec& y,
                                   const arma::ivec& cc,   // 0=obs,1=cens
                                   const int lag,
                                   const arma::mat& adj_matcom, // Ws (nReg x nReg)
                                   const arma::vec& theta,       // [beta(1:p), sigma2, rho_s, rho_t, rho_st]
                                   const arma::vec& lower,       // length = n_cen
                                   const arma::vec& upper,       // length = n_cen
                                   const double jitter = 1e-10,
                                   const double min_prob = 1e-300) {
  // --- unpack theta & basic size checks ---
  const int p = x.n_cols;
  if ((int)theta.n_elem < p + 4) stop("theta must have length >= p+4.");
  arma::vec beta = theta.subvec(0, p-1);
  const double sigma2 = theta(p);
  const double rho_s  = theta(p+1);
  const double rho_t  = theta(p+2);
  const double rho_st = theta(p+3);

  const int nReg = adj_matcom.n_rows;
  if (adj_matcom.n_cols != (unsigned)nReg) stop("adj_matcom must be square.");
  if (lag <= 0) stop("lag must be >= 1.");
  const int nn = nReg * lag;

  if ((int)y.n_elem != nn)
    stop("length(y)=%d but expected nReg*lag=%d*%d=%d.",
         (int)y.n_elem, nReg, lag, nn);

  if ((int)x.n_rows != nn)
    stop("nrow(x)=%d but must equal length(y)=%d.", (int)x.n_rows, (int)y.n_elem);

  if ((int)cc.n_elem != nn)
    stop("length(cc)=%d but must equal length(y)=%d.", (int)cc.n_elem, (int)y.n_elem);

  // ensure cc is only 0 or 1
  if (!find((cc != 0) && (cc != 1)).is_empty())
    stop("cc must contain only 0 (observed) or 1 (censored).");

  // --- indices ---
  arma::uvec idx_obs = find(cc == 0);
  arma::uvec idx_cen = find(cc == 1);
  const std::size_t n_obs = idx_obs.n_elem;
  const std::size_t n_cen = idx_cen.n_elem;

  if ((int)lower.n_elem != (int)n_cen || (int)upper.n_elem != (int)n_cen)
    stop("lower/upper length must equal number of censored entries: n_cen=%d.", (int)n_cen);

  // --- linear predictor ---
  arma::vec eta = x * beta;

  // --- precision Q (no full covariance) ---
  arma::mat Wt = make_Wt(lag);
  arma::mat Q  = buildQST_cpp(adj_matcom, Wt, rho_s, rho_t, rho_st);

  if ((int)Q.n_rows != nn || (int)Q.n_cols != nn)
    stop("Q is %dx%d but expected %dx%d (nReg*lag). Check adj_matcom and lag.",
         (int)Q.n_rows, (int)Q.n_cols, nn, nn);

  // --- partition Q ---
  arma::mat Q_oo, Q_oc, Q_co, Q_cc;
  if (n_obs > 0) Q_oo = Q.submat(idx_obs, idx_obs);
  if (n_obs > 0 && n_cen > 0) {
    Q_oc = Q.submat(idx_obs, idx_cen);
    Q_co = Q.submat(idx_cen, idx_obs);
  }
  if (n_cen > 0) Q_cc = Q.submat(idx_cen, idx_cen);

  // --- Cholesky of Q_cc ---
  arma::mat Lcc;
  if (n_cen > 0) {
    bool okcc = chol(Lcc, Q_cc, "lower");
    if (!okcc && jitter > 0) {
      Q_cc.diag() += jitter;
      okcc = chol(Lcc, Q_cc, "lower");
    }
    if (!okcc) stop("Cholesky failed for Q_cc.");
  }

  // --- B = Q_cc^{-1} Q_co (if needed) ---
  arma::mat B; // (n_cen x n_obs)
  if (n_cen > 0 && n_obs > 0) {
    arma::mat T = solve(trimatl(Lcc), Q_co);    // Lcc * T = Q_co
    B = solve(trimatu(Lcc.t()), T);             // Lcc^T * B = T
  }

  // --- M = Q_oo - Q_oc * B (Schur complement) ---
  arma::mat M;
  if (n_obs > 0) {
    M = Q_oo;
    if (n_cen > 0) M -= Q_oc * B;
  }

  // --- Cholesky of M for observed part ---
  arma::mat Lm;
  if (n_obs > 0) {
    bool okm = chol(Lm, M, "lower");
    if (!okm && jitter > 0) {
      M.diag() += jitter;
      okm = chol(Lm, M, "lower");
    }
    if (!okm) stop("Cholesky failed for M (observed Schur complement).");
  }

  // --- observed loglik ---
  double loglik_obs = 0.0;
  if (n_obs > 0) {
    arma::vec r = y.elem(idx_obs) - eta.elem(idx_obs);
    arma::vec v = Q_oo * r;
    if (n_cen > 0) v -= Q_oc * (B * r);
    const double quad = (1.0 / sigma2) * dot(r, v);
    const double logdetM = 2.0 * sum(log(Lm.diag()));
    const double logdetS = logdetM - (double)n_obs * std::log(sigma2);
    loglik_obs = -0.5 * (double)n_obs * std::log(2.0 * M_PI)
      + 0.5 * logdetS
    - 0.5 * quad;
  }

  // --- censored pmvnorm (pass plain base vectors/matrix) ---
  double loglik_cen = 0.0;
  if (n_cen > 0) {
    arma::vec rc  = (n_obs > 0) ? (y.elem(idx_obs) - eta.elem(idx_obs)) : arma::vec();
    arma::vec muc = eta.elem(idx_cen) - ((n_obs > 0) ? (B * rc) : arma::vec(n_cen, fill::zeros));

    // Sigma_c = sigma2 * Q_cc^{-1}   via Lcc
    arma::mat Icc = eye<mat>(n_cen, n_cen);
    arma::mat X = solve(trimatl(Lcc), Icc);
    arma::mat Qcc_inv = solve(trimatu(Lcc.t()), X);
    arma::mat Sigmac = sigma2 * Qcc_inv;

    // Build plain base R vectors/matrix (no attributes) for mvtnorm
    Environment mvtnorm = Environment::namespace_env("mvtnorm");
    Function pmvnorm = mvtnorm["pmvnorm"];

    const int d = (int)muc.n_elem;
    NumericVector lowerR(d), upperR(d), meanR(d);
    for (int i = 0; i < d; ++i) {
      lowerR[i] = (i < (int)lower.n_elem) ? lower[i] : R_NegInf;
      upperR[i] = (i < (int)upper.n_elem) ? upper[i] : R_PosInf;
      meanR[i]  = muc[i];
    }
    NumericMatrix sigmaR(d, d);
    for (int j = 0; j < d; ++j)
      for (int i = 0; i < d; ++i)
        sigmaR(i, j) = Sigmac(i, j);

    NumericVector probSE = pmvnorm(
      _["lower"] = lowerR,
      _["upper"] = upperR,
      _["mean"]  = meanR,
      _["sigma"] = sigmaR
    // add tolerances if needed: _["abseps"]=1e-8, _["releps"]=1e-6, _["maxpts"]=1e7
    );

    double prob = std::max(min_prob, (double)probSE[0]);
    loglik_cen = std::log(prob);
  }

  return loglik_obs + loglik_cen;
}







// One MVN sample with robust Cholesky (adds tiny jitter only if needed)
inline arma::vec rmvnorm1_chol(const arma::vec& mu, arma::mat Sigma) {
  const std::size_t n = mu.n_rows;
  Sigma = 0.5 * (Sigma + Sigma.t());          // enforce symmetry
  arma::mat L;
  bool ok = arma::chol(L, Sigma, "lower");
  if (!ok) {
    // tiny diagonal jitter; escalate a few times if needed
    for (int k = 0; k < 4 && !ok; ++k) {
      Sigma.diag() += std::pow(10.0, -10 + 2*k);  // 1e-10, 1e-8, 1e-6, 1e-4
      ok = arma::chol(L, Sigma, "lower");
    }
    if (!ok) stop("Cholesky failed (Sigma not SPD).");
  }
  arma::vec z = arma::randn<arma::vec>(n);
  return mu + L * z;
}

//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
// [[Rcpp::export]]
arma::vec pred_dist_spatempcens_ind_cpp(const arma::vec& theta,
                                             Rcpp::Nullable<Rcpp::IntegerVector> predind,
                                             const arma::mat& Xmat,
                                             arma::vec ytotobs,
                                             const int lags,
                                             const arma::ivec& indicens,
                                             arma::mat adjmatinf) {
  const int n = static_cast<int>(Xmat.n_rows);
  const int p = static_cast<int>(Xmat.n_cols);

  if ((int)theta.n_elem < p + 4) {
    stop("theta too short: need at least p+4.");
  }

  // Unpack theta (same layout as en R)
  const arma::vec beta  = theta.subvec(0, p-1);
  const double   sigma2 = theta[p];
  const double   rho    = theta[p+1];
  const double   psi    = theta[p+2];
  const double   phi    = theta[p+3];

  // Insert censored y values
  const arma::uvec idx_cens = arma::find(indicens == 1);
  const int n_cens = static_cast<int>(idx_cens.n_elem);
  if ((int)theta.n_elem < p + 4 + n_cens)
    stop("theta lacks ycensF values to fill censored entries.");

  if (n_cens > 0) {
    const arma::vec ycensF = theta.subvec(p + 4, p + 4 + n_cens - 1);
    ytotobs.elem(idx_cens) = ycensF;
  }

  // Force strict lower-triangular form like in R (upper.tri=0)
  if (adjmatinf.n_rows != adjmatinf.n_cols)
    stop("adjmatinf must be square.");
  for (arma::uword j = 1; j < adjmatinf.n_cols; ++j) {
    for (arma::uword i = 0; i < j; ++i) {
      adjmatinf(i, j) = 0.0;
    }
  }

  // Build covariance (scaled by sigma2)
  arma::mat varcov = sigma2 * spatimecovar_2_cpp(lags, adjmatinf, rho, psi, phi);
  if ((int)varcov.n_rows != n || (int)varcov.n_cols != n)
    stop("spatimecovar_2_cpp returned wrong dimension.");

  // Fast branch: predind == NULL -> sample N(X beta, I)
  if (predind.isNull()) {
    const arma::vec mu = Xmat * beta;
    // I_n covariance => just add standard normal z
    return mu + arma::randn<arma::vec>(n);
  }

  // With predind: conditional Gaussian
  Rcpp::IntegerVector predind_iv(predind.get());
  if (predind_iv.size() != n) stop("predind length must equal nrow(Xmat).");

  const arma::ivec predmask = Rcpp::as<arma::ivec>(predind_iv);
  const arma::uvec idx_obs  = arma::find(predmask == 0);
  const arma::uvec idx_pred = arma::find(predmask == 1);

  if (idx_obs.n_elem == 0u)  stop("No observed indices (predind==0).");
  if (idx_pred.n_elem == 0u) stop("No prediction indices (predind==1).");

  // Means
  const arma::vec mu_full = Xmat * beta;
  const arma::vec mu_obs  = mu_full.elem(idx_obs);
  const arma::vec mu_pred = mu_full.elem(idx_pred);

  // Blocks of Sigma
  const arma::mat Sig_oo = varcov.submat(idx_obs,  idx_obs);
  const arma::mat Sig_pp = varcov.submat(idx_pred, idx_pred);
  const arma::mat Sig_po = varcov.submat(idx_pred, idx_obs);
  const arma::mat Sig_op = varcov.submat(idx_obs,  idx_pred);

  // Cholesky of Sigma_oo
  arma::mat Loo;
  bool ok = arma::chol(Loo, Sig_oo, "lower");
  if (!ok) {
    // light jitter only on the observed block
    arma::mat Sig_oo_j = Sig_oo;
    for (int k = 0; k < 4 && !ok; ++k) {
      Sig_oo_j.diag() += std::pow(10.0, -12 + 2*k); // 1e-12..1e-6
      ok = arma::chol(Loo, Sig_oo_j, "lower");
    }
    if (!ok) stop("Cholesky(Sig_oo) failed.");
  }

  // Compute w = Sigma_oo^{-1} (y_obs - mu_obs) via triangular solves
  arma::vec b  = ytotobs.elem(idx_obs) - mu_obs;
  arma::vec v  = arma::solve(arma::trimatl(Loo), b,  arma::solve_opts::fast);
  arma::vec w  = arma::solve(arma::trimatu(Loo.t()), v, arma::solve_opts::fast);

  // mu_cond = mu_pred + Sig_po * w
  arma::vec mu_cond = mu_pred + Sig_po * w;

  // For Schur complement: X = Sigma_oo^{-1} * Sig_op
  // First Y = L^{-1} * Sig_op, then X = L^{-T} * Y
  arma::mat Y = arma::solve(arma::trimatl(Loo), Sig_op, arma::solve_opts::fast);
  arma::mat X = arma::solve(arma::trimatu(Loo.t()), Y,   arma::solve_opts::fast);

  // Sig_cond = Sig_pp - Sig_po * X
  arma::mat Sig_cond = Sig_pp - Sig_po * X;

  // One conditional sample
  return rmvnorm1_chol(mu_cond, Sig_cond);
}




 // ========================================================================
 // ========== AR(2) version with auto-expansion of short 'indicens' =======
 // ========================================================================
 //' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::vec pred_dist_spatempcens_ar2_ind_cpp(
     const arma::vec& theta,                                 // [beta(1:p), sigma2, rho, psi, phi1, phi2, ycensF...]
     Rcpp::Nullable<Rcpp::IntegerVector> predind,            // length n (0=obs, 1=pred) or NULL
     const arma::mat& Xmat,                                  // n x p
     arma::vec ytotobs,                                      // length n (FULL vector; obs at obs slots)
     const int lags,
     arma::ivec indicens,                                    // length n OR length sum(predind==0) (auto-expanded)
     arma::mat adjmatinf) {                                  // square

   const int n = static_cast<int>(Xmat.n_rows);
   const int p = static_cast<int>(Xmat.n_cols);

   // Basic size checks that do not depend on predind yet
   if ((int)ytotobs.n_elem != n)
     stop("ytotobs length (%d) must equal nrow(Xmat) (%d).", (int)ytotobs.n_elem, n);
   if (adjmatinf.n_rows != adjmatinf.n_cols)
     stop("adjmatinf must be square.");
   if ((int)theta.n_elem < p + 5)
     stop("theta too short: need at least p+5 (beta, sigma2, rho, psi, phi1, phi2).");

   // Unpack theta (same order as your R function)
   const arma::vec beta  = theta.subvec(0, p-1);
   const double   sigma2 = theta[p];
   const double   rho    = theta[p+1];
   const double   psi    = theta[p+2];
   const double   phi1   = theta[p+3];
   const double   phi2   = theta[p+4];

   // If predind is provided, prepare masks and auto-expand a short 'indicens'
   arma::ivec predmask;      // will hold 0/1 split if predind given
   arma::uvec idx_obs, idx_pred;

   if (predind.isNull()) {
     // No predind -> indicens MUST be full length (to place ycensF)
     if ((int)indicens.n_elem != n)
       stop("When predind is NULL, indicens length must equal nrow(Xmat) (%d). Got %d.",
            n, (int)indicens.n_elem);
   } else {
     Rcpp::IntegerVector piv(predind.get());
     if (piv.size() != n)
       stop("predind length (%d) must equal nrow(Xmat) (%d).", piv.size(), n);

     predmask = Rcpp::as<arma::ivec>(piv);
     idx_obs  = arma::find(predmask == 0);
     idx_pred = arma::find(predmask == 1);

     if (idx_obs.n_elem == 0u)  stop("No observed indices (predind==0).");
     if (idx_pred.n_elem == 0u) stop("No prediction indices (predind==1).");

     // --- Auto-expand short indicens (length == n_obs) to full length n ---
     if ((int)indicens.n_elem == (int)idx_obs.n_elem) {
       arma::ivec indic_full(n, arma::fill::zeros);
       indic_full.elem(idx_obs) = indicens;
       indicens = std::move(indic_full);
     } else if ((int)indicens.n_elem != n) {
       stop("indicens length (%d) must equal nrow(Xmat) (%d) or sum(predind==0) (%d).",
            (int)indicens.n_elem, n, (int)idx_obs.n_elem);
     }
   }

   // Insert censored values from tail of theta into ytotobs
   const arma::uvec idx_cens = arma::find(indicens == 1);
   const int n_cens = (int)idx_cens.n_elem;
   if ((int)theta.n_elem < p + 5 + n_cens)
     stop("theta lacks ycensF values: has %d tail, needs %d.",
          (int)theta.n_elem - (p + 5), n_cens);
   if (n_cens > 0) {
     const arma::vec ycensF = theta.subvec(p + 5, p + 5 + n_cens - 1);
     ytotobs.elem(idx_cens) = ycensF;
   }

   // Zero upper triangle of adjmatinf (matches your earlier DAGAR/CAR prep)
   for (arma::uword j = 1; j < adjmatinf.n_cols; ++j)
     for (arma::uword i = 0; i < j; ++i) adjmatinf(i, j) = 0.0;

   // Covariance (scaled)
   arma::mat varcov = sigma2 * spatimecovar_ar2_cpp(lags, adjmatinf, rho, psi, phi1, phi2);
   if ((int)varcov.n_rows != n || (int)varcov.n_cols != n)
     stop("spatimecovar_ar2_cpp returned %dx%d, expected %dx%d.",
          (int)varcov.n_rows, (int)varcov.n_cols, n, n);

   // If no predind: draw from N(X beta, I)
   if (predind.isNull()) {
     const arma::vec mu = Xmat * beta;
     return mu + arma::randn<arma::vec>(n);
   }

   // Conditional Gaussian
   const arma::vec mu_full = Xmat * beta;
   const arma::vec mu_obs  = mu_full.elem(idx_obs);
   const arma::vec mu_pred = mu_full.elem(idx_pred);

   const arma::mat Sig_oo = varcov.submat(idx_obs,  idx_obs);
   const arma::mat Sig_pp = varcov.submat(idx_pred, idx_pred);
   const arma::mat Sig_po = varcov.submat(idx_pred, idx_obs);
   const arma::mat Sig_op = varcov.submat(idx_obs,  idx_pred);

   arma::mat Loo;
   bool ok = arma::chol(Loo, Sig_oo, "lower");
   if (!ok) {
     arma::mat S = Sig_oo;
     for (int k = 0; k < 4 && !ok; ++k) {
       S.diag() += std::pow(10.0, -12 + 2*k);  // 1e-12..1e-6
       ok = arma::chol(Loo, S, "lower");
     }
     if (!ok) stop("Cholesky(Sigma_oo) failed.");
   }

   arma::vec b = ytotobs.elem(idx_obs) - mu_obs;
   arma::vec v = arma::solve(arma::trimatl(Loo), b,  arma::solve_opts::fast);
   arma::vec w = arma::solve(arma::trimatu(Loo.t()), v, arma::solve_opts::fast);

   arma::vec mu_cond = mu_pred + Sig_po * w;

   arma::mat Y = arma::solve(arma::trimatl(Loo), Sig_op, arma::solve_opts::fast);
   arma::mat X = arma::solve(arma::trimatu(Loo.t()), Y,   arma::solve_opts::fast);
   arma::mat Sig_cond = Sig_pp - Sig_po * X;

   return rmvnorm1_chol(mu_cond, Sig_cond);
 }




//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::vec pred_dist_spatempcens_car_ind_cpp(
     const arma::vec& theta,                         // [beta(1:p), sigma2, rho_s, rho_t, rho_st, ycensF...]
     Rcpp::Nullable<Rcpp::IntegerVector> predind,    // length n (0=obs,1=pred) or NULL
     const arma::mat& Xmat,                          // n x p (full)
     arma::vec ytotobs,                              // length n (FULL vector)
     const int lags,
     const arma::ivec& indicens,                     // length n (FULL 0/1 mask)
     const arma::mat& adjmatinf) {                   // spatial adjacency (SYMMETRIC)

   const int n = static_cast<int>(Xmat.n_rows);
   const int p = static_cast<int>(Xmat.n_cols);

   if ((int)theta.n_elem < p + 5)
     stop("theta too short: need at least p+5.");
   if ((int)ytotobs.n_elem != n)
     stop("ytotobs length must equal nrow(Xmat).");
   if ((int)indicens.n_elem != n)
     stop("indicens length must equal nrow(Xmat).");
   if (adjmatinf.n_rows != adjmatinf.n_cols)
     stop("adjmatinf must be square.");

   // ---- unpack theta (CAR layout) ----
   const arma::vec beta  = theta.subvec(0, p-1);
   const double   sigma2 = theta[p];
   const double   rho_s  = theta[p+1];
   const double   rho_t  = theta[p+2];
   const double   rho_st = theta[p+3];

   // insert censored y’s
   const arma::uvec idx_cens = arma::find(indicens == 1);
   const int n_cens = static_cast<int>(idx_cens.n_elem);
   if ((int)theta.n_elem < p + 4 + n_cens)
     stop("theta lacks ycensF values to fill censored entries (need %d).", n_cens);
   if (n_cens > 0) {
     const arma::vec ycensF = theta.subvec(p + 4, p + 4 + n_cens - 1);
     ytotobs.elem(idx_cens) = ycensF;
   }

   // temporal adjacency: |i-j| == 1
   arma::mat Wt(lags, lags, arma::fill::zeros);
   for (int i = 0; i < lags; ++i) {
     if (i+1 < lags) Wt(i, i+1) = 1.0;
     if (i-1 >= 0)   Wt(i, i-1) = 1.0;
   }

   // covariance
   arma::mat K = spatimecovarcar_zero_cpp(adjmatinf, Wt, rho_s, rho_t, rho_st);
   if ((int)K.n_rows != n || (int)K.n_cols != n)
     stop("spatimecovarcar_zero_cpp returned wrong dimension.");
   arma::mat varcov = sigma2 * 0.5 * (K + K.t());   // enforce symmetry

   // no predind: sample N(Xβ, I)
   if (predind.isNull()) {
     const arma::vec mu = Xmat * beta;
     return mu + arma::randn<arma::vec>(n);
   }

   // conditional Gaussian
   Rcpp::IntegerVector predind_iv(predind.get());
   if (predind_iv.size() != n) stop("predind length must equal nrow(Xmat).");
   const arma::ivec predmask = Rcpp::as<arma::ivec>(predind_iv);
   const arma::uvec idx_obs  = arma::find(predmask == 0);
   const arma::uvec idx_pred = arma::find(predmask == 1);
   if (idx_obs.n_elem == 0u)  stop("No observed indices (predind==0).");
   if (idx_pred.n_elem == 0u) stop("No prediction indices (predind==1).");

   const arma::vec mu_full = Xmat * beta;
   const arma::vec mu_obs  = mu_full.elem(idx_obs);
   const arma::vec mu_pred = mu_full.elem(idx_pred);

   arma::mat Sig_oo = varcov.submat(idx_obs,  idx_obs);
   arma::mat Sig_pp = varcov.submat(idx_pred, idx_pred);
   arma::mat Sig_po = varcov.submat(idx_pred, idx_obs);
   arma::mat Sig_op = varcov.submat(idx_obs,  idx_pred);

   // enforce symmetry on the observed block (and jitter if needed)
   Sig_oo = 0.5 * (Sig_oo + Sig_oo.t());
   arma::mat Loo;
   bool ok = arma::chol(Loo, Sig_oo, "lower");
   if (!ok) {
     for (int k = 0; k < 4 && !ok; ++k) {
       Sig_oo.diag() += std::pow(10.0, -12 + 2*k);  // 1e-12..1e-6
       ok = arma::chol(Loo, Sig_oo, "lower");
     }
     if (!ok) stop("Cholesky(Sigma_oo) failed.");
   }

   arma::vec b = ytotobs.elem(idx_obs) - mu_obs;
   arma::vec v = arma::solve(arma::trimatl(Loo), b,  arma::solve_opts::fast);
   arma::vec w = arma::solve(arma::trimatu(Loo.t()), v, arma::solve_opts::fast);

   arma::vec mu_cond = mu_pred + Sig_po * w;

   arma::mat Y = arma::solve(arma::trimatl(Loo), Sig_op, arma::solve_opts::fast);
   arma::mat X = arma::solve(arma::trimatu(Loo.t()), Y,   arma::solve_opts::fast);
   arma::mat Sig_cond = Sig_pp - Sig_po * X;
   Sig_cond = 0.5 * (Sig_cond + Sig_cond.t());      // keep symmetric before sampling

   return rmvnorm1_chol(mu_cond, Sig_cond);
 }



//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 double pred_dens_spatempcens_ind_cpp(const arma::vec& theta,
                                      Rcpp::Nullable<Rcpp::IntegerVector> predind,
                                      const arma::mat& Xmat,
                                      arma::vec ytotobs,
                                      const int lags,
                                      const arma::ivec& indicens,
                                      arma::mat adjmatinf,
                                      const arma::vec& ynew,   // values whose density you want
                                      const bool logd = true)  // return log-density (default)
 {
   const int n = static_cast<int>(Xmat.n_rows);
   const int p = static_cast<int>(Xmat.n_cols);

   if ((int)theta.n_elem < p + 4) {
     Rcpp::stop("theta too short: need at least p+4.");
   }

   // Unpack theta (same layout as in R)
   const arma::vec beta  = theta.subvec(0, p-1);
   const double   sigma2 = theta[p];
   const double   rho    = theta[p+1];
   const double   psi    = theta[p+2];
   const double   phi    = theta[p+3];

   // Fill censored y values from theta
   const arma::uvec idx_cens = arma::find(indicens == 1);
   const int n_cens = static_cast<int>(idx_cens.n_elem);
   if ((int)theta.n_elem < p + 4 + n_cens)
     Rcpp::stop("theta lacks ycensF values to fill censored entries.");

   if (n_cens > 0) {
     const arma::vec ycensF = theta.subvec(p + 4, p + 4 + n_cens - 1);
     ytotobs.elem(idx_cens) = ycensF;
   }

   // Force strict lower-triangular form like in R (upper.tri=0)
   if (adjmatinf.n_rows != adjmatinf.n_cols)
     Rcpp::stop("adjmatinf must be square.");
   for (arma::uword j = 1; j < adjmatinf.n_cols; ++j) {
     for (arma::uword i = 0; i < j; ++i) {
       adjmatinf(i, j) = 0.0;
     }
   }

   // Full covariance (scaled by sigma2)
   arma::mat varcov = sigma2 * spatimecovar_2_cpp(lags, adjmatinf, rho, psi, phi);
   if ((int)varcov.n_rows != n || (int)varcov.n_cols != n)
     Rcpp::stop("spatimecovar_2_cpp returned wrong dimension.");

   const arma::vec mu_full = Xmat * beta;

   // ===== Case A: predind == NULL -> evaluate density under N(mu, I) for all entries =====
   if (predind.isNull()) {
     if ((int)ynew.n_elem != n)
       Rcpp::stop("ynew length must equal nrow(Xmat) when predind is NULL.");

     // N(mu_full, I_n)
     arma::vec r = ynew - mu_full;
     // log |I| = 0; quadratic form is r^T r
     double logdens = -0.5 * ( (double)n * std::log(2.0 * M_PI) + arma::dot(r, r) );
     return logd ? logdens : std::exp(logdens);
   }

   // ===== Case B: conditional density p(y_pred | y_obs) with predind mask =====
   Rcpp::IntegerVector predind_iv(predind.get());
   if (predind_iv.size() != n) Rcpp::stop("predind length must equal nrow(Xmat).");

   const arma::ivec predmask = Rcpp::as<arma::ivec>(predind_iv);
   const arma::uvec idx_obs  = arma::find(predmask == 0);
   const arma::uvec idx_pred = arma::find(predmask == 1);

   if (idx_obs.n_elem == 0u)  Rcpp::stop("No observed indices (predind==0).");
   if (idx_pred.n_elem == 0u) Rcpp::stop("No prediction indices (predind==1).");

   const int k = static_cast<int>(idx_pred.n_elem);
   if ((int)ynew.n_elem != k)
     Rcpp::stop("ynew length must equal number of pred indices (sum(predind==1)).");

   // Partitioned means
   const arma::vec mu_obs  = mu_full.elem(idx_obs);
   const arma::vec mu_pred = mu_full.elem(idx_pred);

   // Partitioned covariances
   const arma::mat Sig_oo = varcov.submat(idx_obs,  idx_obs);
   const arma::mat Sig_pp = varcov.submat(idx_pred, idx_pred);
   const arma::mat Sig_po = varcov.submat(idx_pred, idx_obs);
   const arma::mat Sig_op = varcov.submat(idx_obs,  idx_pred);

   // Cholesky of Sigma_oo with light jitter fallback
   arma::mat Loo;
   bool ok = arma::chol(Loo, Sig_oo, "lower");
   if (!ok) {
     arma::mat Sig_oo_j = Sig_oo;
     for (int t = 0; t < 4 && !ok; ++t) {
       Sig_oo_j.diag() += std::pow(10.0, -12 + 2*t); // 1e-12, 1e-10, 1e-8, 1e-6
       ok = arma::chol(Loo, Sig_oo_j, "lower");
     }
     if (!ok) Rcpp::stop("Cholesky(Sig_oo) failed.");
   }

   // w = Sigma_oo^{-1} (y_obs - mu_obs)
   arma::vec b  = ytotobs.elem(idx_obs) - mu_obs;
   arma::vec v  = arma::solve(arma::trimatl(Loo), b,  arma::solve_opts::fast);
   arma::vec w  = arma::solve(arma::trimatu(Loo.t()), v, arma::solve_opts::fast);

   // Conditional mean: mu_pred + Sig_po * w
   arma::vec mu_cond = mu_pred + Sig_po * w;

   // Schur complement for conditional covariance
   arma::mat Y = arma::solve(arma::trimatl(Loo), Sig_op, arma::solve_opts::fast);
   arma::mat X = arma::solve(arma::trimatu(Loo.t()), Y,   arma::solve_opts::fast);
   arma::mat Sig_cond = Sig_pp - Sig_po * X;

   // Cholesky of Sigma_cond with jitter fallback
   arma::mat Lc;
   ok = arma::chol(Lc, Sig_cond, "lower");
   if (!ok) {
     arma::mat Sig_cj = Sig_cond;
     for (int t = 0; t < 6 && !ok; ++t) {
       Sig_cj.diag() += std::pow(10.0, -12 + 2*t); // up to 1e-2 if really ill-conditioned
       ok = arma::chol(Lc, Sig_cj, "lower");
     }
     if (!ok) Rcpp::stop("Cholesky(Sigma_cond) failed.");
   }

   // Evaluate multivariate normal log-density:
   // log p(y | mu, Sigma) = -0.5*k*log(2π) - sum(log(diag(L))) - 0.5 * || L^{-1}(y-mu) ||^2
   arma::vec r = ynew - mu_cond;
   arma::vec z = arma::solve(arma::trimatl(Lc), r, arma::solve_opts::fast);
   const double quad = arma::dot(z, z);
   const double logdet = arma::sum(arma::log(Lc.diag()));
   const double logdens = -0.5 * ( (double)k * std::log(2.0 * M_PI) + 2.0*logdet + quad );

   return logd ? logdens : std::exp(logdens);
 }



//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 double pred_dens_spatempcens_car_ind_cpp(
     const arma::vec& theta,                         // [beta(1:p), sigma2, rho_s, rho_t, rho_st, ycensF...]
     Rcpp::Nullable<Rcpp::IntegerVector> predind,    // length n (0=obs,1=pred) or NULL
     const arma::mat& Xmat,                          // n x p (full)
     arma::vec ytotobs,                              // length n (FULL vector)
     const int lags,
     const arma::ivec& indicens,                     // length n (FULL 0/1 mask)
     const arma::mat& adjmatinf,                     // spatial adjacency (SYMMETRIC)
     const arma::vec& ynew,                          // values whose density you want
     const bool logd = true)                         // return log-density (default)
 {
   const int n = static_cast<int>(Xmat.n_rows);
   const int p = static_cast<int>(Xmat.n_cols);

   if ((int)theta.n_elem < p + 5)
     Rcpp::stop("theta too short: need at least p+5.");
   if ((int)ytotobs.n_elem != n)
     Rcpp::stop("ytotobs length must equal nrow(Xmat).");
   if ((int)indicens.n_elem != n)
     Rcpp::stop("indicens length must equal nrow(Xmat).");
   if (adjmatinf.n_rows != adjmatinf.n_cols)
     Rcpp::stop("adjmatinf must be square.");

   // ---- unpack theta (CAR layout) ----
   const arma::vec beta  = theta.subvec(0, p-1);
   const double   sigma2 = theta[p];
   const double   rho_s  = theta[p+1];
   const double   rho_t  = theta[p+2];
   const double   rho_st = theta[p+3];

   // insert censored y’s from theta (if any)
   const arma::uvec idx_cens = arma::find(indicens == 1);
   const int n_cens = static_cast<int>(idx_cens.n_elem);
   if ((int)theta.n_elem < p + 4 + n_cens)
     Rcpp::stop("theta lacks ycensF values to fill censored entries (need %d).", n_cens);
   if (n_cens > 0) {
     const arma::vec ycensF = theta.subvec(p + 4, p + 4 + n_cens - 1);
     ytotobs.elem(idx_cens) = ycensF;
   }

   // temporal adjacency: |i-j| == 1
   arma::mat Wt(lags, lags, arma::fill::zeros);
   for (int i = 0; i < lags; ++i) {
     if (i+1 < lags) Wt(i, i+1) = 1.0;
     if (i-1 >= 0)   Wt(i, i-1) = 1.0;
   }

   // Build CAR-based covariance K, then symmetrize & scale
   arma::mat K = spatimecovarcar_zero_cpp(adjmatinf, Wt, rho_s, rho_t, rho_st);
   if ((int)K.n_rows != n || (int)K.n_cols != n)
     Rcpp::stop("spatimecovarcar_zero_cpp returned wrong dimension.");
   arma::mat varcov = sigma2 * 0.5 * (K + K.t());   // enforce symmetry

   const arma::vec mu_full = Xmat * beta;

   // ===== Case A: predind == NULL -> evaluate joint density p(ynew | theta) under N(mu_full, varcov) =====
   if (predind.isNull()) {
     if ((int)ynew.n_elem != n)
       Rcpp::stop("When predind is NULL, ynew length must equal nrow(Xmat).");

     // Cholesky of full covariance with gentle jitter fallback
     arma::mat L;
     bool ok = arma::chol(L, 0.5 * (varcov + varcov.t()), "lower");
     if (!ok) {
       arma::mat Vj = varcov;
       for (int k = 0; k < 6 && !ok; ++k) {
         Vj.diag() += std::pow(10.0, -12 + 2*k); // up to 1e-2 if necessary
         ok = arma::chol(L, 0.5*(Vj + Vj.t()), "lower");
       }
       if (!ok) Rcpp::stop("Cholesky(varcov) failed.");
     }

     arma::vec r = ynew - mu_full;
     arma::vec z = arma::solve(arma::trimatl(L), r, arma::solve_opts::fast);
     const double quad   = arma::dot(z, z);
     const double logdet = arma::sum(arma::log(L.diag()));
     const double logdens = -0.5 * ( (double)n * std::log(2.0 * M_PI) + 2.0*logdet + quad );
     return logd ? logdens : std::exp(logdens);
   }

   // ===== Case B: conditional density p(y_pred | y_obs, theta) with predind mask =====
   Rcpp::IntegerVector predind_iv(predind.get());
   if (predind_iv.size() != n) Rcpp::stop("predind length must equal nrow(Xmat).");

   const arma::ivec predmask = Rcpp::as<arma::ivec>(predind_iv);
   const arma::uvec idx_obs  = arma::find(predmask == 0);
   const arma::uvec idx_pred = arma::find(predmask == 1);
   if (idx_obs.n_elem == 0u)  Rcpp::stop("No observed indices (predind==0).");
   if (idx_pred.n_elem == 0u) Rcpp::stop("No prediction indices (predind==1).");

   const int k = static_cast<int>(idx_pred.n_elem);
   if ((int)ynew.n_elem != k)
     Rcpp::stop("ynew length must equal number of pred indices (sum(predind==1)).");

   // Partitioned means
   const arma::vec mu_obs  = mu_full.elem(idx_obs);
   const arma::vec mu_pred = mu_full.elem(idx_pred);

   // Partitioned covariances
   arma::mat Sig_oo = varcov.submat(idx_obs,  idx_obs);
   arma::mat Sig_pp = varcov.submat(idx_pred, idx_pred);
   arma::mat Sig_po = varcov.submat(idx_pred, idx_obs);
   arma::mat Sig_op = varcov.submat(idx_obs,  idx_pred);

   // Cholesky of Sigma_oo (symmetrize & jitter if needed)
   Sig_oo = 0.5 * (Sig_oo + Sig_oo.t());
   arma::mat Loo;
   bool ok = arma::chol(Loo, Sig_oo, "lower");
   if (!ok) {
     arma::mat Soo = Sig_oo;
     for (int t = 0; t < 6 && !ok; ++t) {
       Soo.diag() += std::pow(10.0, -12 + 2*t); // up to 1e-2 if necessary
       ok = arma::chol(Loo, Soo, "lower");
     }
     if (!ok) Rcpp::stop("Cholesky(Sigma_oo) failed.");
   }

   // w = Sigma_oo^{-1} (y_obs - mu_obs)
   arma::vec b = ytotobs.elem(idx_obs) - mu_obs;
   arma::vec v = arma::solve(arma::trimatl(Loo), b,  arma::solve_opts::fast);
   arma::vec w = arma::solve(arma::trimatu(Loo.t()), v, arma::solve_opts::fast);

   // Conditional mean
   arma::vec mu_cond = mu_pred + Sig_po * w;

   // Conditional covariance via Schur complement
   arma::mat Y = arma::solve(arma::trimatl(Loo), Sig_op, arma::solve_opts::fast);
   arma::mat X = arma::solve(arma::trimatu(Loo.t()), Y,   arma::solve_opts::fast);
   arma::mat Sig_cond = Sig_pp - Sig_po * X;
   Sig_cond = 0.5 * (Sig_cond + Sig_cond.t());

   // Cholesky of Sigma_cond with jitter fallback
   arma::mat Lc;
   ok = arma::chol(Lc, Sig_cond, "lower");
   if (!ok) {
     arma::mat Sc = Sig_cond;
     for (int t = 0; t < 6 && !ok; ++t) {
       Sc.diag() += std::pow(10.0, -12 + 2*t);
       ok = arma::chol(Lc, Sc, "lower");
     }
     if (!ok) Rcpp::stop("Cholesky(Sigma_cond) failed.");
   }

   // Multivariate normal log-density:
   // log p(y | mu, Sigma) = -0.5*k*log(2π) - sum(log(diag(L))) - 0.5 * || L^{-1}(y-mu) ||^2
   arma::vec r = ynew - mu_cond;
   arma::vec z = arma::solve(arma::trimatl(Lc), r, arma::solve_opts::fast);
   const double quad   = arma::dot(z, z);
   const double logdet = arma::sum(arma::log(Lc.diag()));
   const double logdens = -0.5 * ( (double)k * std::log(2.0 * M_PI) + 2.0*logdet + quad );

   return logd ? logdens : std::exp(logdens);
 }



//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 double pred_dens_spatempcens_ar2_ind_cpp(
     const arma::vec& theta,                                 // [beta(1:p), sigma2, rho, psi, phi1, phi2, ycensF...]
     Rcpp::Nullable<Rcpp::IntegerVector> predind,            // length n (0=obs, 1=pred) or NULL
     const arma::mat& Xmat,                                  // n x p
     arma::vec ytotobs,                                      // length n (FULL vector; obs at obs slots)
     const int lags,
     arma::ivec indicens,                                    // length n OR length sum(predind==0) (auto-expanded)
     arma::mat adjmatinf,                                    // square
     const arma::vec& ynew,                                  // values whose density to evaluate
     const bool logd = true)                                 // return log-density (default)
 {
   const int n = static_cast<int>(Xmat.n_rows);
   const int p = static_cast<int>(Xmat.n_cols);

   // Basic checks
   if ((int)ytotobs.n_elem != n)
     Rcpp::stop("ytotobs length (%d) must equal nrow(Xmat) (%d).", (int)ytotobs.n_elem, n);
   if (adjmatinf.n_rows != adjmatinf.n_cols)
     Rcpp::stop("adjmatinf must be square.");
   if ((int)theta.n_elem < p + 5)
     Rcpp::stop("theta too short: need at least p+5 (beta, sigma2, rho, psi, phi1, phi2).");

   // Unpack theta
   const arma::vec beta  = theta.subvec(0, p-1);
   const double   sigma2 = theta[p];
   const double   rho    = theta[p+1];
   const double   psi    = theta[p+2];
   const double   phi1   = theta[p+3];
   const double   phi2   = theta[p+4];

   // predind handling + auto-expand short indicens
   arma::ivec predmask; arma::uvec idx_obs, idx_pred;
   bool has_pred = !predind.isNull();

   if (!has_pred) {
     if ((int)indicens.n_elem != n)
       Rcpp::stop("When predind is NULL, indicens length must equal nrow(Xmat) (%d). Got %d.",
                  n, (int)indicens.n_elem);
   } else {
     Rcpp::IntegerVector piv(predind.get());
     if (piv.size() != n)
       Rcpp::stop("predind length (%d) must equal nrow(Xmat) (%d).", piv.size(), n);

     predmask = Rcpp::as<arma::ivec>(piv);
     idx_obs  = arma::find(predmask == 0);
     idx_pred = arma::find(predmask == 1);
     if (idx_obs.n_elem == 0u)  Rcpp::stop("No observed indices (predind==0).");
     if (idx_pred.n_elem == 0u) Rcpp::stop("No prediction indices (predind==1).");

     if ((int)indicens.n_elem == (int)idx_obs.n_elem) {
       arma::ivec indic_full(n, arma::fill::zeros);
       indic_full.elem(idx_obs) = indicens;
       indicens = std::move(indic_full);
     } else if ((int)indicens.n_elem != n) {
       Rcpp::stop("indicens length (%d) must equal nrow(Xmat) (%d) or sum(predind==0) (%d).",
                  (int)indicens.n_elem, n, (int)idx_obs.n_elem);
     }
   }

   // Insert censored values from tail of theta into ytotobs
   const arma::uvec idx_cens = arma::find(indicens == 1);
   const int n_cens = (int)idx_cens.n_elem;
   if ((int)theta.n_elem < p + 5 + n_cens)
     Rcpp::stop("theta lacks ycensF values: has %d tail, needs %d.",
                (int)theta.n_elem - (p + 5), n_cens);
   if (n_cens > 0) {
     const arma::vec ycensF = theta.subvec(p + 5, p + 5 + n_cens - 1);
     ytotobs.elem(idx_cens) = ycensF;
   }

   // Zero upper triangle of adjmatinf (as in your pipeline)
   for (arma::uword j = 1; j < adjmatinf.n_cols; ++j)
     for (arma::uword i = 0; i < j; ++i) adjmatinf(i, j) = 0.0;

   // Covariance (scaled)
   arma::mat varcov = sigma2 * spatimecovar_ar2_cpp(lags, adjmatinf, rho, psi, phi1, phi2);
   if ((int)varcov.n_rows != n || (int)varcov.n_cols != n)
     Rcpp::stop("spatimecovar_ar2_cpp returned %dx%d, expected %dx%d.",
                (int)varcov.n_rows, (int)varcov.n_cols, n, n);
   varcov = 0.5 * (varcov + varcov.t()); // enforce symmetry

   const arma::vec mu_full = Xmat * beta;

   // =================== Case A: No mask -> joint density p(ynew | theta) ===================
   if (!has_pred) {
     if ((int)ynew.n_elem != n)
       Rcpp::stop("When predind is NULL, ynew length must equal nrow(Xmat) (%d).", n);

     arma::mat L;
     bool ok = arma::chol(L, varcov, "lower");
     if (!ok) {
       arma::mat Vj = varcov;
       for (int k = 0; k < 6 && !ok; ++k) {
         Vj.diag() += std::pow(10.0, -12 + 2*k); // 1e-12 ... 1e-2
         ok = arma::chol(L, 0.5*(Vj + Vj.t()), "lower");
       }
       if (!ok) Rcpp::stop("Cholesky(varcov) failed.");
     }

     arma::vec r = ynew - mu_full;
     arma::vec z = arma::solve(arma::trimatl(L), r, arma::solve_opts::fast);
     const double quad   = arma::dot(z, z);
     const double logdet = arma::sum(arma::log(L.diag()));
     const double logdens = -0.5 * ( (double)n * std::log(2.0 * M_PI) + 2.0*logdet + quad );
     return logd ? logdens : std::exp(logdens);
   }

   // ===== Case B: Conditional density p(y_pred | y_obs, theta) with predind mask =====
   const int k = static_cast<int>(idx_pred.n_elem);
   if ((int)ynew.n_elem != k)
     Rcpp::stop("ynew length (%d) must equal number of pred indices (sum(predind==1)=%d).",
                (int)ynew.n_elem, k);

   // Partitioned means and covariances
   const arma::vec mu_obs  = mu_full.elem(idx_obs);
   const arma::vec mu_pred = mu_full.elem(idx_pred);

   arma::mat Sig_oo = varcov.submat(idx_obs,  idx_obs);
   arma::mat Sig_pp = varcov.submat(idx_pred, idx_pred);
   arma::mat Sig_po = varcov.submat(idx_pred, idx_obs);
   arma::mat Sig_op = varcov.submat(idx_obs,  idx_pred);

   // Cholesky of Sigma_oo (symmetrize & jitter if needed)
   Sig_oo = 0.5 * (Sig_oo + Sig_oo.t());
   arma::mat Loo;
   bool ok = arma::chol(Loo, Sig_oo, "lower");
   if (!ok) {
     arma::mat Soo = Sig_oo;
     for (int t = 0; t < 6 && !ok; ++t) {
       Soo.diag() += std::pow(10.0, -12 + 2*t);
       ok = arma::chol(Loo, Soo, "lower");
     }
     if (!ok) Rcpp::stop("Cholesky(Sigma_oo) failed.");
   }

   // w = Sigma_oo^{-1} (y_obs - mu_obs)
   arma::vec b = ytotobs.elem(idx_obs) - mu_obs;
   arma::vec v = arma::solve(arma::trimatl(Loo), b,  arma::solve_opts::fast);
   arma::vec w = arma::solve(arma::trimatu(Loo.t()), v, arma::solve_opts::fast);

   // Conditional mean
   arma::vec mu_cond = mu_pred + Sig_po * w;

   // Conditional covariance via Schur complement
   arma::mat Y = arma::solve(arma::trimatl(Loo), Sig_op, arma::solve_opts::fast);
   arma::mat X = arma::solve(arma::trimatu(Loo.t()), Y,   arma::solve_opts::fast);
   arma::mat Sig_cond = Sig_pp - Sig_po * X;
   Sig_cond = 0.5 * (Sig_cond + Sig_cond.t());

   // Cholesky of Sigma_cond with jitter fallback
   arma::mat Lc;
   ok = arma::chol(Lc, Sig_cond, "lower");
   if (!ok) {
     arma::mat Sc = Sig_cond;
     for (int t = 0; t < 6 && !ok; ++t) {
       Sc.diag() += std::pow(10.0, -12 + 2*t);
       ok = arma::chol(Lc, 0.5*(Sc + Sc.t()), "lower");
     }
     if (!ok) Rcpp::stop("Cholesky(Sigma_cond) failed.");
   }

   // Multivariate normal log-density:
   // log p(y | mu, Sigma) = -0.5*k*log(2π) - sum(log(diag(L))) - 0.5 * || L^{-1}(y-mu) ||^2
   arma::vec r = ynew - mu_cond;
   arma::vec z = arma::solve(arma::trimatl(Lc), r, arma::solve_opts::fast);
   const double quad   = arma::dot(z, z);
   const double logdet = arma::sum(arma::log(Lc.diag()));
   const double logdens = -0.5 * ( (double)k * std::log(2.0 * M_PI) + 2.0*logdet + quad );

   return logd ? logdens : std::exp(logdens);
 }


//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::vec pred_pointwise_logdens_ar2_ind_cpp(
     const arma::vec& theta,                  // [beta(1:p), sigma2, rho, psi, phi1, phi2, ycensF...]
     const Rcpp::IntegerVector& predind,      // length n (0=obs, 1=pred)  -- REQUIRED
     const arma::mat& Xmat,                   // n x p
     arma::vec ytotobs,                       // length n (FULL vector; obs at obs slots)
     const int lags,
     arma::ivec indicens,                     // length n OR length sum(predind==0) (auto-expanded)
     arma::mat adjmatinf,                     // square
     const arma::vec& ynew                    // length k (values at pred positions)
 ){
   const int n = (int)Xmat.n_rows, p = (int)Xmat.n_cols;
   if ((int)predind.size() != n) Rcpp::stop("predind length must equal n");
   if (adjmatinf.n_rows != adjmatinf.n_cols) Rcpp::stop("adjmatinf must be square");
   if ((int)theta.n_elem < p + 5) Rcpp::stop("theta too short");

   // Unpack theta
   const arma::vec beta  = theta.subvec(0, p-1);
   const double   sigma2 = theta[p];
   const double   rho    = theta[p+1];
   const double   psi    = theta[p+2];
   const double   phi1   = theta[p+3];
   const double   phi2   = theta[p+4];

   // Masks
   arma::ivec pm = Rcpp::as<arma::ivec>(predind);
   arma::uvec idx_obs  = arma::find(pm == 0);
   arma::uvec idx_pred = arma::find(pm == 1);
   if (idx_obs.n_elem == 0u || idx_pred.n_elem == 0u) Rcpp::stop("Need both obs and pred indices");

   const int k = (int)idx_pred.n_elem;
   if ((int)ynew.n_elem != k) Rcpp::stop("ynew length must equal sum(predind==1)");

   // Auto-expand short indicens if provided only for observed part
   if ((int)indicens.n_elem == (int)idx_obs.n_elem) {
     arma::ivec indic_full(n, arma::fill::zeros);
     indic_full.elem(idx_obs) = indicens;
     indicens = std::move(indic_full);
   } else if ((int)indicens.n_elem != n) {
     Rcpp::stop("indicens must have length n or sum(predind==0)");
   }

   // Fill censored y's from theta tail
   const arma::uvec idx_cens = arma::find(indicens == 1);
   const int n_cens = (int)idx_cens.n_elem;
   if ((int)theta.n_elem < p + 5 + n_cens)
     Rcpp::stop("theta lacks ycensF values");
   if (n_cens > 0) {
     const arma::vec ycensF = theta.subvec(p + 5, p + 5 + n_cens - 1);
     ytotobs.elem(idx_cens) = ycensF;
   }

   // Zero upper triangle (match your pipeline)
   for (arma::uword j = 1; j < adjmatinf.n_cols; ++j)
     for (arma::uword i = 0; i < j; ++i) adjmatinf(i, j) = 0.0;

   // Covariance & mean
   arma::mat varcov = sigma2 * spatimecovar_ar2_cpp(lags, adjmatinf, rho, psi, phi1, phi2);
   if ((int)varcov.n_rows != n || (int)varcov.n_cols != n)
     Rcpp::stop("spatimecovar_ar2_cpp returned wrong size");
   varcov = 0.5*(varcov + varcov.t()); // symmetrize for safety

   const arma::vec mu_full = Xmat * beta;
   const arma::vec mu_obs  = mu_full.elem(idx_obs);
   const arma::vec mu_pred = mu_full.elem(idx_pred);

   // Blocks
   arma::mat Sig_oo = varcov.submat(idx_obs,  idx_obs);
   arma::mat Sig_pp = varcov.submat(idx_pred, idx_pred);
   arma::mat Sig_po = varcov.submat(idx_pred, idx_obs);
   arma::mat Sig_op = varcov.submat(idx_obs,  idx_pred);
   Sig_oo = 0.5*(Sig_oo + Sig_oo.t());

   // Cholesky of Sig_oo
   arma::mat Loo;
   bool ok = arma::chol(Loo, Sig_oo, "lower");
   if (!ok) {
     arma::mat S = Sig_oo;
     for (int t=0; t<6 && !ok; ++t) { S.diag() += std::pow(10.0, -12 + 2*t); ok = arma::chol(Loo, S, "lower"); }
     if (!ok) Rcpp::stop("Cholesky(Sig_oo) failed");
   }

   // w = Sig_oo^{-1}(y_obs - mu_obs)
   arma::vec b = ytotobs.elem(idx_obs) - mu_obs;
   arma::vec v = arma::solve(arma::trimatl(Loo), b,  arma::solve_opts::fast);
   arma::vec w = arma::solve(arma::trimatu(Loo.t()), v, arma::solve_opts::fast);

   // Conditional mean
   arma::vec mu_cond = mu_pred + Sig_po * w;

   // For diag(Sig_cond):
   // Sig_cond = Sig_pp - Sig_po * (Sig_oo^{-1} Sig_op)
   arma::mat Y = arma::solve(arma::trimatl(Loo), Sig_op, arma::solve_opts::fast);   // Loo^{-1} Sig_op
   arma::mat X = arma::solve(arma::trimatu(Loo.t()), Y,   arma::solve_opts::fast);   // Sig_oo^{-1} Sig_op
   // diag(Sig_po * X) = rowwise sum of (Sig_po % X^T)
   arma::vec diag_term = arma::sum(Sig_po % X.t(), 1);
   arma::vec var_cond  = Sig_pp.diag() - diag_term;

   // Safety clamp (numerical jitter if tiny negatives due to roundoff)
   const double eps = 1e-12;
   for (int i=0; i<k; ++i) if (var_cond(i) < eps) var_cond(i) = eps;

   // Univariate normal log-densities (pointwise)
   arma::vec z = ynew - mu_cond;
   arma::vec logdens = -0.5 * ( arma::log(2.0 * M_PI * var_cond) + (z % z) / var_cond );
   return logdens;
 }







//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::vec pred_pointwise_logdens_ar1_ind_cpp(
     const arma::vec& theta,                  // [beta(1:p), sigma2, rho, psi, phi, ycensF...]
     const Rcpp::IntegerVector& predind,      // length n (0=obs, 1=pred)  -- REQUIRED
     const arma::mat& Xmat,                   // n x p
     arma::vec ytotobs,                       // length n (FULL vector; obs at obs slots)
     const int lags,
     arma::ivec indicens,                     // length n OR length sum(predind==0) (auto-expanded)
     arma::mat adjmatinf,                     // square
     const arma::vec& ynew                    // length k (values at pred positions)
 ){
   const int n = (int)Xmat.n_rows, p = (int)Xmat.n_cols;
   if ((int)predind.size() != n) Rcpp::stop("predind length must equal n");
   if (adjmatinf.n_rows != adjmatinf.n_cols) Rcpp::stop("adjmatinf must be square");
   if ((int)theta.n_elem < p + 4) Rcpp::stop("theta too short (need beta, sigma2, rho, psi, phi)");

   // Unpack theta (AR1)
   const arma::vec beta  = theta.subvec(0, p-1);
   const double   sigma2 = theta[p];
   const double   rho    = theta[p+1];
   const double   psi    = theta[p+2];
   const double   phi    = theta[p+3];

   // Masks
   arma::ivec pm = Rcpp::as<arma::ivec>(predind);
   arma::uvec idx_obs  = arma::find(pm == 0);
   arma::uvec idx_pred = arma::find(pm == 1);
   if (idx_obs.n_elem == 0u || idx_pred.n_elem == 0u) Rcpp::stop("Need both obs and pred indices");

   const int k = (int)idx_pred.n_elem;
   if ((int)ynew.n_elem != k) Rcpp::stop("ynew length must equal sum(predind==1)");

   // Auto-expand short indicens if given only for observed part
   if ((int)indicens.n_elem == (int)idx_obs.n_elem) {
     arma::ivec indic_full(n, arma::fill::zeros);
     indic_full.elem(idx_obs) = indicens;
     indicens = std::move(indic_full);
   } else if ((int)indicens.n_elem != n) {
     Rcpp::stop("indicens must have length n or sum(predind==0)");
   }

   // Fill censored y's from theta tail
   const arma::uvec idx_cens = arma::find(indicens == 1);
   const int n_cens = (int)idx_cens.n_elem;
   if ((int)theta.n_elem < p + 4 + n_cens)
     Rcpp::stop("theta lacks ycensF values");
   if (n_cens > 0) {
     const arma::vec ycensF = theta.subvec(p + 4, p + 4 + n_cens - 1);
     ytotobs.elem(idx_cens) = ycensF;
   }

   // Zero upper triangle (your convention)
   for (arma::uword j = 1; j < adjmatinf.n_cols; ++j)
     for (arma::uword i = 0; i < j; ++i) adjmatinf(i, j) = 0.0;

   // Covariance & mean  (swap builder name if needed)
   arma::mat varcov = sigma2 * spatimecovar_2_cpp(lags, adjmatinf, rho, psi, phi);
   // If your code uses spatimecovar_2_cpp(...), use that instead:
   // arma::mat varcov = sigma2 * spatimecovar_2_cpp(lags, adjmatinf, rho, psi, phi);

   if ((int)varcov.n_rows != n || (int)varcov.n_cols != n)
     Rcpp::stop("AR1 varcov wrong size");
   varcov = 0.5*(varcov + varcov.t()); // symmetrize

   const arma::vec mu_full = Xmat * beta;
   const arma::vec mu_obs  = mu_full.elem(idx_obs);
   const arma::vec mu_pred = mu_full.elem(idx_pred);

   // Blocks
   arma::mat Sig_oo = varcov.submat(idx_obs,  idx_obs);
   arma::mat Sig_pp = varcov.submat(idx_pred, idx_pred);
   arma::mat Sig_po = varcov.submat(idx_pred, idx_obs);
   arma::mat Sig_op = varcov.submat(idx_obs,  idx_pred);
   Sig_oo = 0.5*(Sig_oo + Sig_oo.t());

   // Cholesky of Sig_oo
   arma::mat Loo;
   bool ok = arma::chol(Loo, Sig_oo, "lower");
   if (!ok) {
     arma::mat S = Sig_oo;
     for (int t=0; t<6 && !ok; ++t) { S.diag() += std::pow(10.0, -12 + 2*t); ok = arma::chol(Loo, S, "lower"); }
     if (!ok) Rcpp::stop("Cholesky(Sig_oo) failed");
   }

   // w = Sig_oo^{-1}(y_obs - mu_obs)
   arma::vec b = ytotobs.elem(idx_obs) - mu_obs;
   arma::vec v = arma::solve(arma::trimatl(Loo), b,  arma::solve_opts::fast);
   arma::vec w = arma::solve(arma::trimatu(Loo.t()), v, arma::solve_opts::fast);

   // Conditional mean
   arma::vec mu_cond = mu_pred + Sig_po * w;

   // diag(Sig_cond) without forming the full k×k matrix
   arma::mat Y = arma::solve(arma::trimatl(Loo), Sig_op, arma::solve_opts::fast);   // Loo^{-1} Sig_op
   arma::mat X = arma::solve(arma::trimatu(Loo.t()), Y,   arma::solve_opts::fast);   // Sig_oo^{-1} Sig_op
   arma::vec diag_term = arma::sum(Sig_po % X.t(), 1);
   arma::vec var_cond  = Sig_pp.diag() - diag_term;

   const double eps = 1e-12;
   for (int i=0; i<k; ++i) if (var_cond(i) < eps) var_cond(i) = eps;

   // Univariate Normal log-densities
   arma::vec z = ynew - mu_cond;
   arma::vec logdens = -0.5 * ( arma::log(2.0 * M_PI * var_cond) + (z % z) / var_cond );
   return logdens;
 }


//' @useDynLib auxiliarcpp, .registration = TRUE
 //' @importFrom Rcpp sourceCpp
 //' @export
 // [[Rcpp::export]]
 arma::vec pred_pointwise_logdens_car_ind_cpp(
     const arma::vec& theta,                  // [beta(1:p), sigma2, rho_s, rho_t, rho_st, ycensF...]
     const Rcpp::IntegerVector& predind,      // length n (0=obs, 1=pred)  -- REQUIRED
     const arma::mat& Xmat,                   // n x p
     arma::vec ytotobs,                       // length n (FULL vector)
     const int lags,
     arma::ivec indicens,                     // length n OR length sum(predind==0) (auto-expanded)
     arma::mat adjmatinf,                     // spatial adjacency (SYMMETRIC)
     const arma::vec& ynew                    // length k (values at pred positions)
 ){
   const int n = (int)Xmat.n_rows, p = (int)Xmat.n_cols;
   if ((int)predind.size() != n) Rcpp::stop("predind length must equal n");
   if (adjmatinf.n_rows != adjmatinf.n_cols) Rcpp::stop("adjmatinf must be square");
   if ((int)theta.n_elem < p + 4) Rcpp::stop("theta too short (need beta, sigma2, rho_s, rho_t, rho_st)");

   // Unpack theta (CAR)
   const arma::vec beta  = theta.subvec(0, p-1);
   const double   sigma2 = theta[p];
   const double   rho_s  = theta[p+1];
   const double   rho_t  = theta[p+2];
   const double   rho_st = theta[p+3];

   // Masks
   arma::ivec pm = Rcpp::as<arma::ivec>(predind);
   arma::uvec idx_obs  = arma::find(pm == 0);
   arma::uvec idx_pred = arma::find(pm == 1);
   if (idx_obs.n_elem == 0u || idx_pred.n_elem == 0u) Rcpp::stop("Need both obs and pred indices");

   const int k = (int)idx_pred.n_elem;
   if ((int)ynew.n_elem != k) Rcpp::stop("ynew length must equal sum(predind==1)");

   // Auto-expand short indicens if given only for observed part
   if ((int)indicens.n_elem == (int)idx_obs.n_elem) {
     arma::ivec indic_full(n, arma::fill::zeros);
     indic_full.elem(idx_obs) = indicens;
     indicens = std::move(indic_full);
   } else if ((int)indicens.n_elem != n) {
     Rcpp::stop("indicens must have length n or sum(predind==0)");
   }

   // Fill censored y's from theta tail
   const arma::uvec idx_cens = arma::find(indicens == 1);
   const int n_cens = (int)idx_cens.n_elem;
   if ((int)theta.n_elem < p + 4 + n_cens)
     Rcpp::stop("theta lacks ycensF values");
   if (n_cens > 0) {
     const arma::vec ycensF = theta.subvec(p + 4, p + 4 + n_cens - 1);
     ytotobs.elem(idx_cens) = ycensF;
   }

   // Temporal adjacency Wt for |i-j|==1
   arma::mat Wt(lags, lags, arma::fill::zeros);
   for (int i = 0; i < lags; ++i) {
     if (i+1 < lags) Wt(i, i+1) = 1.0;
     if (i-1 >= 0)   Wt(i, i-1) = 1.0;
   }

   // Build CAR covariance
   arma::mat K = spatimecovarcar_zero_cpp(adjmatinf, Wt, rho_s, rho_t, rho_st);
   if ((int)K.n_rows != n || (int)K.n_cols != n)
     Rcpp::stop("spatimecovarcar_zero_cpp returned wrong size");
   arma::mat varcov = sigma2 * 0.5 * (K + K.t()); // enforce symmetry

   const arma::vec mu_full = Xmat * beta;
   const arma::vec mu_obs  = mu_full.elem(idx_obs);
   const arma::vec mu_pred = mu_full.elem(idx_pred);

   // Blocks
   arma::mat Sig_oo = varcov.submat(idx_obs,  idx_obs);
   arma::mat Sig_pp = varcov.submat(idx_pred, idx_pred);
   arma::mat Sig_po = varcov.submat(idx_pred, idx_obs);
   arma::mat Sig_op = varcov.submat(idx_obs,  idx_pred);
   Sig_oo = 0.5*(Sig_oo + Sig_oo.t());

   // Cholesky of Sig_oo
   arma::mat Loo;
   bool ok = arma::chol(Loo, Sig_oo, "lower");
   if (!ok) {
     arma::mat S = Sig_oo;
     for (int t=0; t<6 && !ok; ++t) { S.diag() += std::pow(10.0, -12 + 2*t); ok = arma::chol(Loo, S, "lower"); }
     if (!ok) Rcpp::stop("Cholesky(Sig_oo) failed");
   }

   // w = Sig_oo^{-1}(y_obs - mu_obs)
   arma::vec b = ytotobs.elem(idx_obs) - mu_obs;
   arma::vec v = arma::solve(arma::trimatl(Loo), b,  arma::solve_opts::fast);
   arma::vec w = arma::solve(arma::trimatu(Loo.t()), v, arma::solve_opts::fast);

   // Conditional mean
   arma::vec mu_cond = mu_pred + Sig_po * w;

   // diag(Sig_cond) without forming full
   arma::mat Y = arma::solve(arma::trimatl(Loo), Sig_op, arma::solve_opts::fast);
   arma::mat X = arma::solve(arma::trimatu(Loo.t()), Y,   arma::solve_opts::fast);
   arma::vec diag_term = arma::sum(Sig_po % X.t(), 1);
   arma::vec var_cond  = Sig_pp.diag() - diag_term;

   const double eps = 1e-12;
   for (int i=0; i<k; ++i) if (var_cond(i) < eps) var_cond(i) = eps;

   // Univariate Normal log-densities
   arma::vec z = ynew - mu_cond;
   arma::vec logdens = -0.5 * ( arma::log(2.0 * M_PI * var_cond) + (z % z) / var_cond );
   return logdens;
 }
