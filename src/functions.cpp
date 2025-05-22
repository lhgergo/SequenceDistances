#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// Function to transcode one amino acid to a corresponding number
int findLetterPosition(char target) {
  std::string aas_unq = "ACDEFGHIKLMNPQRSTVWY";
  std::vector<char> vec;
  vec.insert(vec.end(), aas_unq.begin(), aas_unq.end());
  int n = vec.size();
  
  for (int p = 0; p < n; ++p) {
    if (vec[p] == target) {
      return p;  // Adding 1 because R indexing starts from 1
    }
  }
  
  // If the letter is not found, return -1
  return -1;
}

// Function to calculate Grantham distance between two amino acids
int grantham_distance(char aa1, char aa2) {
  // Grantham matrix values for amino acid properties
  int grantham_matrix[20][20] = {
    {0, 195, 126, 107, 113, 60, 86, 94, 106, 96, 84, 111, 27, 91, 112, 99, 58, 64, 148, 112},
    {195, 0, 154, 170, 205, 159, 174, 198, 202, 198, 196, 139, 169, 154, 180, 112, 149, 192, 215, 194},
    {126, 154, 0, 45, 177, 94, 81, 168, 101, 172, 160, 23, 108, 61, 96, 65, 85, 152, 181, 160},
    {107, 170, 45, 0, 140, 98, 40, 134, 56, 138, 126, 42, 93, 29, 54, 80, 65, 121, 152, 122},
    {113, 205, 177, 140, 0, 153, 100, 21, 102, 22, 28, 158, 114, 116, 97, 155, 103, 50, 40, 22},
    {60, 159, 94, 98, 153, 0, 98, 135, 127, 138, 127, 80, 42, 87, 125, 56, 59, 109, 184, 147},
    {86, 174, 81, 40, 100, 98, 0, 94, 32, 99, 87, 68, 77, 24, 29, 89, 47, 84, 115, 83},
    {94, 198, 168, 134, 21, 135, 94, 0, 102, 5, 10, 149, 95, 109, 97, 142, 89, 29, 61, 33},
    {106, 202, 101, 56, 102, 127, 32, 102, 0, 107, 95, 94, 103, 53, 26, 121, 78, 97, 110, 85},
    {96, 198, 172, 138, 22, 138, 99, 5, 107, 0, 15, 153, 98, 113, 102, 145, 92, 32, 61, 36},
    {84, 196, 160, 126, 28, 127, 87, 10, 95, 15, 0, 142, 87, 101, 91, 135, 81, 21, 67, 36},
    {111, 139, 23, 42, 158, 80, 68, 149, 94, 153, 142, 0, 91, 46, 86, 46, 65, 133, 174, 143},
    {27, 169, 108, 93, 114, 42, 77, 95, 103, 98, 87, 91, 0, 76, 103, 74, 38, 68, 147, 110},
    {91, 154, 61, 29, 116, 87, 24, 109, 53, 113, 101, 46, 76, 0, 43, 68, 42, 96, 130, 99},
    {112, 180, 96, 54, 97, 125, 29, 97, 26, 102, 91, 86, 103, 43, 0, 110, 71, 96, 101, 77},
    {99, 112, 65, 80, 155, 56, 89, 142, 121, 145, 135, 46, 74, 68, 110, 0, 58, 124, 177, 144},
    {58, 149, 85, 65, 103, 59, 47, 89, 78, 92, 81, 65, 38, 42, 71, 58, 0, 69, 128, 92},
    {64, 192, 152, 121, 50, 109, 84, 29, 97, 32, 21, 133, 68, 96, 96, 124, 69, 0, 88, 55},
    {148, 215, 181, 152, 40, 184, 115, 61, 110, 61, 67, 174, 147, 130, 101, 177, 128, 88, 0, 37},
    {112, 194, 160, 122, 22, 147, 83, 33, 85, 36, 36, 143, 110, 99, 77, 144, 92, 55, 37, 0}
  };
  
  // Convert amino acids to index
  int index1 = findLetterPosition(aa1);
  int index2 = findLetterPosition(aa2);
  
  // Return the Grantham distance
  return grantham_matrix[index1][index2];
}

// Function to calculate Grantham distances between two sets of peptides
// [[Rcpp::export]]
NumericMatrix GranthamDistance(StringVector peptides1, StringVector peptides2) {
  int n1 = peptides1.size();
  int n2 = peptides2.size();
  
  NumericMatrix distances(n1, n2);
  
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      // Get peptide sequences as strings
      std::vector<char> seq1_vec;
      std::vector<char> seq2_vec;
      seq1_vec.insert(seq1_vec.end(), peptides1[i].begin(), peptides1[i].end());
      seq2_vec.insert(seq2_vec.end(), peptides2[j].begin(), peptides2[j].end());
      
      // Calculate Grantham distance for each pair of amino acids in the sequences
      int distance = 0;
      for (int k = 0; k < seq1_vec.size() && k < seq2_vec.size(); ++k) {
        distance += grantham_distance(seq1_vec[k], seq2_vec[k]);
      }
      
      // Store the calculated distance
      distances(i, j) = distance;
    }
  }
  
  return distances;
}

#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate BLOSUM62 similarity between two amino acids
int blosum62_similarity(char aa1, char aa2) {
  int blosum62[20][20] = {{4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2},
                          {0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2},
                          {-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3},
                          {-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2},
                          {-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3},
                          {0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3},
                          {-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2},
                          {-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1},
                          {-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2},
                          {-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1},
                          {-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1},
                          {-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2},
                          {-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3},
                          {-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1},
                          {-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2},
                          {1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2},
                          {0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2},
                          {0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1},
                          {-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2},
                          {-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7}};
  int pos1 = findLetterPosition(aa1);
  int pos2 = findLetterPosition(aa2);
  if (pos1 == -1 || pos2 == -1) return -999;
  return blosum62[pos1][pos2];
}


// Function to calculate BLOSUM62 similarities between two sets of peptides
// [[Rcpp::export]]
NumericMatrix BLOSUM62sim(StringVector peptides1, StringVector peptides2) {
  int n1 = peptides1.size();
  int n2 = peptides2.size();
  
  NumericMatrix similarities(n1, n2);
  
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      // Get peptide sequences as strings
      std::vector<char> seq1_vec;
      std::vector<char> seq2_vec;
      seq1_vec.insert(seq1_vec.end(), peptides1[i].begin(), peptides1[i].end());
      seq2_vec.insert(seq2_vec.end(), peptides2[j].begin(), peptides2[j].end());
      
      // Calculate Grantham similarity for each pair of amino acids in the sequences
      int similarity = 0;
      for (int k = 0; k < seq1_vec.size() && k < seq2_vec.size(); ++k) {
        similarity += blosum62_similarity(seq1_vec[k], seq2_vec[k]);
      }
      
      // Store the calculated similarity
      similarities(i, j) = similarity;
    }
  }
  
  return similarities;
}