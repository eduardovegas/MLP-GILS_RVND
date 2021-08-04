#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <chrono>
#include "readData.h"

#define IMAX 10
#define MAXDB 99999999999.99

using namespace std;

double ** matrizAdj; // matriz de adjacencia
int dimension; // quantidade total de vertices

template <typename A, typename B>
void zip(const vector<A> &a, const vector<B> &b, vector<pair<A,B> > &zipped)
{
  int size = a.size();
  for(size_t i = 0; i < size; ++i)
  {
    zipped.push_back(make_pair(a[i], b[i]));
  }
}

template <typename A, typename B>
void unzip(const vector<pair<A, B> > &zipped, vector<A> &a, vector<B> &b)
{
  int size = a.size();
  for(size_t i = 0; i < size; i++)
  {
    a[i] = zipped[i].first;
    b[i] = zipped[i].second;
  }
}

void printData()
{
  cout << "dimension: " << dimension << endl;
  for (int i = 1; i <= dimension; i++)
  {
    for (int j = 1; j <= dimension; j++)
    {
      cout << matrizAdj[i][j] << " ";
    }
    cout << endl;
  }
}

void printReopt(vector<vector<int>>& W, vector<vector<double>>& T, vector<vector<double>>& C, int selector)
{
  if(selector == 1)
  {
    cout << "W : " << endl;
    for(int i = 0; i < W.size(); i++)
    {
      for(int j = 0; j < W[i].size(); j++)
      {
        cout << W[i][j] << " ";
      }
      cout << endl;
    }
  };

  if(selector == 2)
  {
    cout << "C : " << endl;
    for(int i = 0; i < C.size(); i++)
    {
      for(int j = 0; j < C[i].size(); j++)
      {
        cout << C[i][j] << " ";
      }
      cout << endl;
    }
  };

  if(selector == 3)
  {
    cout << "T : " << endl;
    for(int i = 0; i < T.size(); i++)
    {
      for(int j = 0; j < T[i].size(); j++)
      {
        cout << T[i][j] << " ";
      }
      cout << endl;
    }
  };
}

void printSolution(vector<int>& sol, double cost, vector<vector<int>>& W, vector<vector<double>>& T, vector<vector<double>>& C)
{
  for(int i = 0; i < sol.size(); i++)
  {
    cout << sol[i] << " ";
  }
  cout << "cost = " << cost << endl;

  cost = 0;
  for(int i = 0; i < sol.size()-1; i++)
  {
    for(int j = 0; j < i; j++)
    {
      cost += matrizAdj[sol[j]][sol[j+1]];
    }
  }
  cout << "calculated cost = " << cost + T[0][dimension] << endl;
  cout << "calculated cost = " << C[0][dimension] + T[0][dimension] << "\n" << endl;
}

void Swap(vector<int>& sol, double* cost, bool* improved, vector<vector<int>>& W, vector<vector<double>>& T, vector<vector<double>>& C)
{
  int c1, t1;
  int c2, t2;
  int c3, t3;
  int c4, t4;
  int best_i;
  int best_j;
  double cheaper_cost = *cost;

  int solSize = sol.size();

  for(int i = 1; i < solSize - 2; i++)
  {
    //primeiro testa swap adjacente ao no i
    int k = i+1;
    
    //a solucao eh dividida nas seguintes subsequencias:
    //(0, 1, ..., i-1) (i, k) (k+1, k+2, ..., 0)

    //subsequencia de (0, 1, ..., i-1) + subsequencia de (k, i) = s1 -> (0, 1, ..., i-1, k, i)
    t1 = T[0][i-1] + matrizAdj[sol[i-1]][sol[k]];
    c1 = C[0][i-1] + 2*(t1) + C[k][i];
    t1 += T[k][i];

    //subsequencia de s1 + subsequencia de (k+1, k+2, ..., 0) = s2 -> (0, 1, ..., i-1, k, i, k+1, ..., 0)
    t2 = t1 + matrizAdj[sol[i]][sol[k+1]];
    //c2 = c1 + (dimension-k-1)*(t2) + C[k+1][dimension];
    c2 = c1 + W[k+1][dimension]*(t2) + C[k+1][dimension]; //W[k+1][dimension] == (dimension-k-1)
    t2 += T[k+1][dimension];

    c2 += t2; //adicionando o T de todo o percurso

    if(c2 < cheaper_cost)
    {
      best_i = i;
      best_j = k;
      *improved = true;
      cheaper_cost = c2;
    }

    for(int j = i + 2; j < solSize - 1; j++) //depois testa todas as outras possibilidades de swap
    {
      //a solucao eh dividida nas seguintes subsequencias:
      //(0, 1, ..., i-1) (i) (i+1, i+2, ..., j-2, j-1) (j) (j+1, j+2, ..., 0)

      //subsequencia de (0, 1, ..., i-1) + subsequencia de (j) = s1 -> (0, 1, ..., i-1, j)
      t1 = T[0][i-1] + matrizAdj[sol[i-1]][sol[j]];
      c1 = C[0][i-1] + 1*(t1) + 0;
      t1 += 0;
      
      //subsequencia de s1 + subsequencia de (i+1, i+2, ..., j-2, j-1) = s2 -> (0, ..., i-1, j, i+1, ..., j-2, j-1)
      t2 = t1 + matrizAdj[sol[j]][sol[i+1]];
      //c2 = c1 + (j-i-1)*(t2) + C[i+1][j-1];
      c2 = c1 + W[i+1][j-1]*(t2) + C[i+1][j-1]; //W[i+1][j-1] == (j-i-1)
      t2 += T[i+1][j-1];
      
      //subsequencia de s2 + subsequencia de (i) = s3 -> (0, ..., i-1, j, i+1, ..., j-2, j-1, i)
      t3 = t2 + matrizAdj[sol[j-1]][sol[i]];
      c3 = c2 + 1*(t3) + 0;
      t3 += 0;

      //subsequencia de s3 + subsequencia de (j+1, j+2, ..., 0) = s4 -> (0, ..., i-1, j, i+1, ..., j-2, j-1, i, j+1, j+2, ..., 0))
      t4 = t3 + matrizAdj[sol[i]][sol[j+1]];
      //c4 = c3 + (dimension-j-1)*(t4) + C[j+1][dimension];
      c4 = c3 + W[j+1][dimension]*(t4) + C[j+1][dimension]; //W[j+1][dimension] == (dimension-j-1)
      t4 += T[j+1][dimension];

      c4 += t4; //adicionando o T de todo o percurso

      if(c4 < cheaper_cost)
      {
        best_i = i;
        best_j = j;
        *improved = true;
        cheaper_cost = c4;
      }
    }
  }

  if(*improved)
  {
    //printf("best_i = %d\nbest_j = %d\n", best_i, best_j);
    *cost = cheaper_cost;
    swap(sol[best_i], sol[best_j]);
  }
}

void opt2(vector<int>& sol, double* cost, bool* improved, vector<vector<int>>& W, vector<vector<double>>& T, vector<vector<double>>& C)
{
  int c1, t1;
  int c2, t2;
  int best_i;
  int best_j;
  double cheaper_cost = *cost;

  int solSize = sol.size();

  for(int i = 1; i < solSize - 4; i++)
  {
    for(int j = i + 3; j < solSize - 1; j++)
    {
      //a solucao eh dividida nas seguintes subsequencias:
      //(0, 1, ..., i-1) (i, i+1, ..., j-1, j) (j+1, j+2, ..., 0)

      //subsequencia de (0, 1, ..., i-1) + subsequencia de (j, j-1, ..., i+1, i) = s1 -> (0, 1, ..., i-1, j, j-1, ..., i+1, i)
      t1 = T[0][i-1] + matrizAdj[sol[i-1]][sol[j]];
      //c1 = C[0][i-1] + (j-i+1)*(t1) + C[j][i];
      c1 = C[0][i-1] + W[j][i]*(t1) + C[j][i]; //W[j][i] == (j-i+1)
      t1 += T[j][i];

      //subsequencia de s1 + subsequencia de (j+1, j+2, ..., 0) = s2 -> (0, 1, ..., i-1, j, j-1, ..., i+1, i, j+1, ..., 0)
      t2 = t1 + matrizAdj[sol[i]][sol[j+1]];
      //c2 = c1 + (dimension-j-1)*(t2) + C[j+1][dimension];
      c2 = c1 + W[j+1][dimension]*(t2) + C[j+1][dimension]; //W[j+1][dimension] == (dimension-j-1)
      t2 += T[j+1][dimension];

      c2 += t2; //adicionando o T de todo o percurso

      if(c2 < cheaper_cost)
      {
        best_i = i;
        best_j = j;
        *improved = true;
        cheaper_cost = c2;
      }
    }
  }

  if(*improved)
  {
    //printf("best_i = %d  best_j = %d\n", best_i, best_j);
    *cost = cheaper_cost;
    reverse(sol.begin() + best_i, sol.begin() + best_j + 1);
  }
}

void reinsertion(vector<int>& sol, double* cost, bool* improved, vector<vector<int>>& W, vector<vector<double>>& T, vector<vector<double>>& C)
{
  int c1, t1;
  int c2, t2;
  int c3, t3;
  int best_i;
  int best_j;
  double cheaper_cost = *cost;

  int solSize = sol.size();

  for(int i = 1; i < solSize - 1; i++)
  {
    for(int j = i+2; j < solSize - 1; j++) //primeiro todas as possibilidades com j > i
    {
      //a solucao eh dividida nas seguintes subsequencias:
      //(0, 1, ..., i-1) (i) (i+1, i+2, ..., j-1, j) (j+1, j+2, ..., 0)

      //subsequencia de (0, 1, ..., i-1) + subsequencia de (i+1, i+2, ..., j-1, j) = s1 -> (0, 1, ..., i-1, i+1, ..., j-1, j)
      t1 = T[0][i-1] + matrizAdj[sol[i-1]][sol[i+1]];
      //c1 = C[0][i-1] + (j-i)*(t1) + C[i+1][j];
      c1 = C[0][i-1] + W[i+1][j]*(t1) + C[i+1][j]; //W[i+1][j] == (j-i)
      t1 += T[i+1][j];

      //subsequencia de s1 + subsequencia de (i) = s2 -> (0, 1, ..., i-1, i+1, ..., j-1, j, i)
      t2 = t1 + matrizAdj[sol[j]][sol[i]];
      c2 = c1 + 1*(t2) + 0;
      t2 += 0;

      //subsequencia de s2 + subsequencia de (j+1, j+2, ..., 0) = s3 -> (0, 1, ..., i-1, i+1, ..., j-1, j, i, j+1, ..., 0)
      t3 = t2 + matrizAdj[sol[i]][sol[j+1]];
      //c3 = c2 + (dimension-j-1)*(t3) + C[j+1][dimension];
      c3 = c2 + W[j+1][dimension]*(t3) + C[j+1][dimension]; //W[j+1][dimension] == (dimension-j-1)
      t3 += T[j+1][dimension];

      c3 += t3; //adicionando o T de todo o percurso

      if(c3 < cheaper_cost)
      {
        best_i = i;
        best_j = j;
        *improved = true;
        cheaper_cost = c3;
      }
    }

    for(int j = 0; j <= i-3; j++) //depois todas as possibilidades com i > j
    {
      //a solucao eh dividida nas seguintes subsequencias:
      //(0, 1, ..., j-1, j) (j+1, j+2, ..., i-2, i-1) (i) (i+1, i+2, ..., 0)

      //subsequencia de (0, 1, ..., j-1, j) + subsequencia de (i) = s1 -> (0, 1, ..., j-1, j, i)
      t1 = T[0][j] + matrizAdj[sol[j]][sol[i]];
      c1 = C[0][j] + 1*(t1) + 0;
      t1 += 0;

      //subsequencia de s1 + subsequencia de (j+1, j+2, ..., i-2, i-1) = s2 -> (0, 1, ..., j-1, j, i, j+1, ..., i-1)
      t2 = t1 + matrizAdj[sol[i]][sol[j+1]];
      //c2 = c1 + (i-j-1)*(t2) + C[j+1][i-1];
      c2 = c1 + W[j+1][i-1]*(t2) + C[j+1][i-1]; //W[j+1][i-1] == (i-j-1)
      t2 += T[j+1][i-1];

      //subsequencia de s2 + subsequencia de (i+1, i+2, ..., 0) = s3 -> (0, 1, ..., j-1, j, i, j+1, ..., i-1, i+1, ..., 0)
      t3 = t2 + matrizAdj[sol[i-1]][sol[i+1]];
      //c3 = c2 + (dimension-i-1)*(t3) + C[i+1][dimension];
      c3 = c2 + W[i+1][dimension]*(t3) + C[i+1][dimension]; //W[i+1][dimension] == (dimension-i-1)
      t3 += T[i+1][dimension];

      c3 += t3; //adicionando o T de todo o percurso

      if(c3 < cheaper_cost)
      {
        best_i = i;
        best_j = j;
        *improved = true;
        cheaper_cost = c3;
      }
    }
  }

  if(*improved)
  {
    //printf("best_i = %d\nbest_j = %d\ndecrease of = %lf\n", sol[best_i], sol[best_j], best_decrease);
    *cost = cheaper_cost;

    int node = sol[best_i];

    if(best_i > best_j)
    {
      sol.erase(sol.begin() + best_i, sol.begin() + best_i + 1);
      sol.insert(sol.begin() + best_j + 1, node);
    }
    else
    {
      sol.erase(sol.begin() + best_i, sol.begin() + best_i + 1);
      sol.insert(sol.begin() + best_j, node);
    }
  }
}

void reinsertion2(vector<int>& sol, double* cost, bool* improved, vector<vector<int>>& W, vector<vector<double>>& T, vector<vector<double>>& C)
{
  int c1, t1;
  int c2, t2;
  int c3, t3;
  int best_i;
  int best_j;
  double cheaper_cost = *cost;

  int solSize = sol.size();

  for(int i = 1; i < solSize - 2; i++)
  {
    for(int j = i+2; j < solSize - 1; j++) //primeiro todas as possibilidades com j > i
    {
      //a solucao eh dividida nas seguintes subsequencias:
      //(0, 1, ..., i-1) (i, i+1) (i+2, i+3, ..., j-1, j) (j+1, j+2, ..., 0)

      //subsequencia de (0, 1, ..., i-1) + subsequencia de (i+2, i+3, ..., j-1, j) = s1 -> (0, 1, ..., i-1, i+2, ..., j-1, j)
      t1 = T[0][i-1] + matrizAdj[sol[i-1]][sol[i+2]];
      //c1 = C[0][i-1] + (j-i-1)*(t1) + C[i+2][j];
      c1 = C[0][i-1] + W[i+2][j]*(t1) + C[i+2][j]; //W[i+2][j] == (j-i-1)
      t1 += T[i+2][j];

      //subsequencia de s1 + subsequencia de (i, i+1) = s2 -> (0, 1, ..., i-1, i+2, ..., j-1, j, i, i+1)
      t2 = t1 + matrizAdj[sol[j]][sol[i]];
      c2 = c1 + 2*(t2) + C[i][i+1];
      t2 += T[i][i+1];

      //subsequencia de s2 + subsequencia de (j+1, j+2, ..., 0) = s3 -> (0, 1, ..., i-1, i+2, ..., j-1, j, i, i+1, j+1, ..., 0)
      t3 = t2 + matrizAdj[sol[i+1]][sol[j+1]];
      //c3 = c2 + (dimension-j-1)*(t3) + C[j+1][dimension];
      c3 = c2 + W[j+1][dimension]*(t3) + C[j+1][dimension]; //W[j+1][dimension] == (dimension-j-1)
      t3 += T[j+1][dimension];

      c3 += t3; //adicionando o T de todo o percurso

      if(c3 < cheaper_cost)
      {
        best_i = i;
        best_j = j;
        *improved = true;
        cheaper_cost = c3;
      }
    }

    for(int j = 0; j <= i-2; j++) //depois todas as possibilidades com i > j
    {
      //a solucao eh dividida nas seguintes subsequencias:
      //(0, 1, ..., j-1, j) (j+1, j+2, ..., i-2, i-1) (i, i+1) (i+2, i+3, ..., 0)

      //subsequencia de (0, 1, ..., j-1, j) + subsequencia de (i, i+1) = s1 -> (0, 1, ..., j-1, j, i, i+1)
      t1 = T[0][j] + matrizAdj[sol[j]][sol[i]];
      c1 = C[0][j] + 2*(t1) + C[i][i+1];
      t1 += T[i][i+1];

      //subsequencia de s1 + subsequencia de (j+1, j+2, ..., i-2, i-1) = s2 -> (0, 1, ..., j-1, j, i, i+1. j+1, ..., i-1)
      t2 = t1 + matrizAdj[sol[i+1]][sol[j+1]];
      //c2 = c1 + (i-j-1)*(t2) + C[j+1][i-1];
      c2 = c1 + W[j+1][i-1]*(t2) + C[j+1][i-1]; //W[j+1][i-1] == (i-j-1)
      t2 += T[j+1][i-1];

      //subsequencia de s2 + subsequencia de (i+2, i+3, ..., 0) = s3 -> (0, 1, ..., j-1, j, i, j+1, ..., i-1, i+2, ..., 0)
      t3 = t2 + matrizAdj[sol[i-1]][sol[i+2]];
      //c3 = c2 + (dimension-i-2)*(t3) + C[i+2][dimension];
      c3 = c2 + W[i+2][dimension]*(t3) + C[i+2][dimension]; //W[i+2][dimension] == (dimension-i-2)
      t3 += T[i+2][dimension];

      c3 += t3; //adicionando o T de todo o percurso

      if(c3 < cheaper_cost)
      {
        best_i = i;
        best_j = j;
        *improved = true;
        cheaper_cost = c3;
      }
    }
  }

  if(*improved)
  {
    //printf("best_i = %d\nbest_j = %d\ndecrease of = %lf\n", sol[best_i], sol[best_j], best_decrease);
    *cost = cheaper_cost;

    vector<int> subsequence(sol.begin() + best_i, sol.begin() + best_i + 2);

    if(best_i > best_j)
    {
      sol.erase(sol.begin() + best_i, sol.begin() + best_i + 2);
      sol.insert(sol.begin() + best_j + 1, subsequence.begin(), subsequence.end());
    }
    else
    {
      sol.erase(sol.begin() + best_i, sol.begin() + best_i + 2);
      sol.insert(sol.begin() + best_j - 1, subsequence.begin(), subsequence.end());
    }
  }
}

void reinsertion3(vector<int>& sol, double* cost, bool* improved, vector<vector<int>>& W, vector<vector<double>>& T, vector<vector<double>>& C)
{
  int c1, t1;
  int c2, t2;
  int c3, t3;
  int best_i;
  int best_j;
  double cheaper_cost = *cost;

  int solSize = sol.size();

  for(int i = 1; i < solSize - 3; i++)
  {
    for(int j = i+3; j < solSize - 1; j++) //primeiro todas as possibilidades com j > i
    {
      //a solucao eh dividida nas seguintes subsequencias:
      //(0, 1, ..., i-1) (i, i+1, i+2) (i+3, i+4, ..., j-1, j) (j+1, j+2, ..., 0)

      //subsequencia de (0, 1, ..., i-1) + subsequencia de (i+3, i+4, ..., j-1, j) = s1 -> (0, 1, ..., i-1, i+3, ..., j-1, j)
      t1 = T[0][i-1] + matrizAdj[sol[i-1]][sol[i+3]];
      //c1 = C[0][i-1] + (j-i-2)*(t1) + C[i+3][j];
      c1 = C[0][i-1] + W[i+3][j]*(t1) + C[i+3][j]; //W[i+3][j] == (j-i-2)
      t1 += T[i+3][j];

      //subsequencia de s1 + subsequencia de (i, i+1, i+2) = s2 -> (0, 1, ..., i-1, i+3, ..., j-1, j, i, i+1, i+2)
      t2 = t1 + matrizAdj[sol[j]][sol[i]];
      c2 = c1 + 3*(t2) + C[i][i+2];
      t2 += T[i][i+2];

      //subsequencia de s2 + subsequencia de (j+1, j+2, ..., 0) = s3 -> (0, 1, ..., i-1, i+3, ..., j-1, j, i, i+1, i+2, j+1, ..., 0)
      t3 = t2 + matrizAdj[sol[i+2]][sol[j+1]];
      //c3 = c2 + (dimension-j-1)*(t3) + C[j+1][dimension];
      c3 = c2 + W[j+1][dimension]*(t3) + C[j+1][dimension]; //W[j+1][dimension] == (dimension-j-1)
      t3 += T[j+1][dimension];

      c3 += t3; //adicionando o T de todo o percurso

      if(c3 < cheaper_cost)
      {
        best_i = i;
        best_j = j;
        *improved = true;
        cheaper_cost = c3;
      }
    }

    for(int j = 0; j <= i-2; j++) //depois todas as possibilidades com i > j
    {
      //a solucao eh dividida nas seguintes subsequencias:
      //(0, 1, ..., j-1, j) (j+1, j+2, ..., i-2, i-1) (i, i+1, i+2) (i+3, i+4, ..., 0)

      //subsequencia de (0, 1, ..., j-1, j) + subsequencia de (i, i+1, i+2) = s1 -> (0, 1, ..., j-1, j, i, i+1, i+2)
      t1 = T[0][j] + matrizAdj[sol[j]][sol[i]];
      c1 = C[0][j] + 3*(t1) + C[i][i+2];
      t1 += T[i][i+2];

      //subsequencia de s1 + subsequencia de (j+1, j+2, ..., i-2, i-1) = s2 -> (0, 1, ..., j-1, j, i, i+1, i+2, j+1, ..., i-1)
      t2 = t1 + matrizAdj[sol[i+2]][sol[j+1]];
      //c2 = c1 + (i-j-1)*(t2) + C[j+1][i-1];
      c2 = c1 + W[j+1][i-1]*(t2) + C[j+1][i-1]; //W[j+1][i-1] == (i-j-1)
      t2 += T[j+1][i-1];

      //subsequencia de s2 + subsequencia de (i+3, i+4, ..., 0) = s3 -> (0, 1, ..., j-1, j, i, i+1, i+2, j+1, ..., i-1, i+3, ..., 0)
      t3 = t2 + matrizAdj[sol[i-1]][sol[i+3]];
      //c3 = c2 + (dimension-i-3)*(t3) + C[i+3][dimension];
      c3 = c2 + W[i+3][dimension]*(t3) + C[i+3][dimension]; //W[i+3][dimension] == (dimension-i-3)
      t3 += T[i+3][dimension];

      c3 += t3; //adicionando o T de todo o percurso

      if(c3 < cheaper_cost)
      {
        best_i = i;
        best_j = j;
        *improved = true;
        cheaper_cost = c3;
      }
    }
  }

  if(*improved)
  {
    //printf("best_i = %d\nbest_j = %d\ndecrease of = %lf\n", sol[best_i], sol[best_j], best_decrease);
    *cost = cheaper_cost;

    vector<int> subsequence(sol.begin() + best_i, sol.begin() + best_i + 3);

    if(best_i > best_j)
    {
      sol.erase(sol.begin() + best_i, sol.begin() + best_i + 3);
      sol.insert(sol.begin() + best_j + 1, subsequence.begin(), subsequence.end());
    }
    else
    {
      sol.erase(sol.begin() + best_i, sol.begin() + best_i + 3);
      sol.insert(sol.begin() + best_j - 2, subsequence.begin(), subsequence.end());
    }
  }
}

void constructW(vector<vector<int>>& W)
{
  int solSize = dimension+1;

  for(int i = 0; i < dimension; i++)
  {
    W[i][i] = 1;
  }
  W[0][0] = 0;
  W[dimension][dimension] = 0;

  //iremos construir o vetor W apenas uma vez no programa, pois ele nao muda
  for(int n = 2; n <= solSize; n++) //numero de cidades que separam as subsequencias, comecando em 2
  {
    int tam = (dimension-n)+2;
    for(int i = 0, j; i < tam; i++) //i eh uma subsequencia de 1 no
    {
      j = (i+n)-1; //j eh outra subsequencia de 1 no que vai concatenar com i

      W[i][j] = W[i][j-1] + W[j][j];

      W[j][i] = W[i][j];
    }
  }
}

void constructReopt(vector<vector<int>>& W, vector<vector<double>>& T, vector<vector<double>>& C, vector<int>& sol)
{
  int solSize = sol.size();
  double tempo_aux = 0;

  for(int i = 0; i < dimension; i++)
  {
    T[i][i] = 0;
    C[i][i] = 0;
  }

  for(int n = 2; n <= solSize; n++) //numero de cidades que separam as subsequencias, comecando em 2
  {
    int tam = (dimension-n)+2;
    for(int i = 0, j; i < tam; i++) //i eh uma subsequencia de 1 no
    {
      j = (i+n)-1; //j eh outra subsequencia de 1 no que vai concatenar com i
      
      tempo_aux = matrizAdj[sol[j-1]][sol[j]]; //por ser uma matriz em que t[i][j] == t[j][i], podemos escolher qualquer um para usar nos dois casos

      T[i][j] = T[i][j-1] + tempo_aux;
      C[i][j] = C[i][j-1] + W[j][j]*(T[i][j]); // + C[j][j] = 0

      T[j][i] = T[i][j];
      C[j][i] = C[j-1][i] + W[j-1][i]*(tempo_aux); // + C[j][j] = 0 e + T[j][j] = 0
    }
  }
}

void construction(vector<int>& sol, double* cost, float alfa)
{
  sol.push_back(1);
  sol.push_back(1);

  vector<int> candidates;
  for(int i = 1; i < dimension; i++)
  {  
    candidates.push_back(i+1);
  }

  int index = 0;
  double custo_aux = 0;

  while(candidates.size() != 0)
  {
    vector<double> deltas;
    vector<int> noInserido;
    vector<int> arestaRemovida;

    for(int i = index, j = i+1; i < sol.size()-1; i++, j++)
    {
      custo_aux = matrizAdj[sol[i]][sol[j]];

      for(int k = 0; k < candidates.size(); k++)
      {
        deltas.push_back(matrizAdj[sol[i]][candidates[k]] + matrizAdj[candidates[k]][sol[j]] - custo_aux);
        noInserido.push_back(k);
        arestaRemovida.push_back(i);
      }

      index = i+1;
    }

    vector<pair<double,int> > zipped;
    vector<pair<double,int> > zipped2;
    zip(deltas, noInserido, zipped);
    zip(deltas, arestaRemovida, zipped2);

    sort(zipped.begin(), zipped.end()); //ordenar os nos inseridos pelos deltas
    sort(zipped2.begin(), zipped2.end()); //ordenar as arestas removidas pelos deltas

    unzip(zipped, deltas, noInserido);
    unzip(zipped2, deltas, arestaRemovida);

    int aux = (int)(floor(alfa*candidates.size()));
    int r = 0;

    if(aux != 0)
    {
      r = rand() % aux;
    }
    
    sol.insert(sol.begin() + (arestaRemovida[r]+1), candidates[noInserido[r]]);
    candidates.erase(candidates.begin() + noInserido[r]);
  }
}

void rvnd(vector<int>& sol, double* cost, vector<vector<int>>& W, vector<vector<double>>& T, vector<vector<double>>& C)
{
  int r = 0;
  bool improved = false;

  vector<int> NL;
  for(int i = 0; i < 5; i++)
  {
    NL.push_back(i);
  }

  while(NL.size() != 0)
  {

    r = rand() % NL.size();

    switch(NL[r])
    {
      case 0:
        Swap(sol, cost, &improved, W, T, C);
        break;
        
      case 1:
        opt2(sol, cost, &improved, W, T, C);
        break;

      case 2:
        reinsertion(sol, cost, &improved, W, T, C);
        break;

      case 3:
        reinsertion2(sol, cost, &improved, W, T, C);
        break;

      case 4:
        reinsertion3(sol, cost, &improved, W, T, C);
        break;
    }

    if(improved)
    {
      NL.clear();
      for(int i = 0; i < 5; i++)
      {
        NL.push_back(i);
      }

      improved = false;

      constructReopt(W, T, C, sol); //atualizar estruturas auxiliares
    }
    else 
    {
      NL.erase(NL.begin() + r);
    }
  }
}

void dbridge(vector<int>& sol, double* cost)
{
  int tam1 = 0; //tam subsequencia 1
  int tam2 = 0; //tam subsequencia 2
  int index1 = 0; //indice de onde comeca subsequencia 1
  int index2 = 0; //indice de onde comeca subsequencia 2

  //tam das subsequencias entre 2 e |V|/10
  if(dimension/10 <= 2)
  {
    tam1 = 2;
    tam2 = 2;
  }
  else
  {
    tam1 = rand() % ((dimension/10)-1) + 2;
    tam2 = rand() % ((dimension/10)-1) + 2;
  }
  
  //indice de inicio da primeira subsequencia dentro de um intervalo valido
  index1 = rand() % ((sol.size()-2)-(tam1-1)) + 1; //rand() % (sol.size()-2) + 1 da um indice dentro vetor excluindo origem e destino

  //procura indice de inicio da segunda subsequencia dentro de um intervalo valido
  while(1)
  {
    index2 = rand() % ((sol.size()-2)-(tam2-1)) + 1; //rand() % (sol.size()-2) + 1 da um indice dentro vetor excluindo origem e destino

    if((index1 < index2) && (index1+tam1 > index2)) //subsequencias com intersecao
      continue;

    if((index2 < index1) && (index2+tam2 > index1)) //subsequencias com intersecao
      continue;

    if(index1 == index2)
      continue;
    
    break;
  }

  if(index2 < index1)
  {
    int aux_index = index1;
    index1 = index2;
    index2 = aux_index;

    int aux_tam = tam1;
    tam1 = tam2;
    tam2 = aux_tam;
  }

  // if(index1 + tam1 == index2)
  // {
  //   *cost = *cost + (matrizAdj[sol[index1-1]][sol[index2]] + matrizAdj[sol[index2+tam2-1]][sol[index1]] + matrizAdj[sol[index1+tam1-1]][sol[index2+tam2]] - matrizAdj[sol[index1-1]][sol[index1]] - matrizAdj[sol[index1+tam1-1]][sol[index1+tam1]] - matrizAdj[sol[index2+tam2-1]][sol[index2+tam2]]);
  // }
  // else
  // {
  //   *cost = *cost + (matrizAdj[sol[index1-1]][sol[index2]] + matrizAdj[sol[index2+tam2-1]][sol[index1+tam1]] + matrizAdj[sol[index2-1]][sol[index1]] + matrizAdj[sol[index1+tam1-1]][sol[index2+tam2]] - matrizAdj[sol[index1-1]][sol[index1]] - matrizAdj[sol[index1+tam1-1]][sol[index1+tam1]] - matrizAdj[sol[index2-1]][sol[index2]] - matrizAdj[sol[index2+tam2-1]][sol[index2+tam2]]);
  // }

  vector<int> subsequence1(sol.begin() + index1, sol.begin() + index1 + tam1);
  vector<int> subsequence2(sol.begin() + index2, sol.begin() + index2 + tam2);

  //printf("best_i = %d\ntam_i = %d\nbest_j = %d\ntam_j = %d\nnew cost = %lf\n", sol[index1], tam1, sol[index2], tam2, *cost);
  
  sol.erase(sol.begin() + index1, sol.begin() + index1 + tam1); //remove a primeira subsequencia
  sol.insert(sol.begin() + index1, subsequence2.begin(), subsequence2.end()); //insere a segunda no lugar

  index2 += tam2 - tam1; //indice da segunda subsequencia deslocado

  sol.erase(sol.begin() + index2, sol.begin() + index2 + tam2); //remove a segunda subsequencia
  sol.insert(sol.begin() + index2, subsequence1.begin(), subsequence1.end()); //insere a primeira no lugar
}

int main(int argc, char** argv)
{
  readData(argc, argv, &dimension, &matrizAdj);

  float alfa;
  int iterIls = 0;
  int Iils = (dimension >= 100) ? 100 : dimension;

  double current_cost;
  double better_cost;
  double best_cost = MAXDB;
  vector<int> current_solution;
  vector<int> better_solution;
  vector<int> best_solution;
  vector<vector<int>> W(dimension+1, vector<int>(dimension+1));
  vector<vector<double>> T(dimension+1, vector<double>(dimension+1));
  vector<vector<double>> C(dimension+1, vector<double>(dimension+1));

  srand(time(NULL));

  //GILS-RVND
  auto timerStart = chrono::system_clock::now();
  
  constructW(W); //a estrutura auxiliar W so precisa ser construida uma vez

  for(int i = 0; i < IMAX; i++)
  {
    iterIls = 0;
    alfa = rand() % 101;
    alfa = alfa/100.0;

    current_solution.clear();

    construction(current_solution, &current_cost, alfa); //solucao inicial
    constructReopt(W, T, C, current_solution); //estruturas auxiliares
    current_cost = C[0][dimension] + T[0][dimension]; //custo inicial
    
    better_cost = current_cost;
    better_solution = current_solution;

    while(iterIls < Iils)
    {
      rvnd(current_solution, &current_cost, W, T, C); //busca por otimo local

      if(current_cost < better_cost)
      {
        iterIls = 0;
        better_cost = current_cost;
        better_solution = current_solution;
      }
      else
      {
        current_solution = better_solution;
        current_cost = better_cost;
      }

      dbridge(current_solution, &current_cost); //pertuba a solucao

      constructReopt(W, T, C, current_solution); //atualiza estruturas auxiliares
      current_cost = C[0][dimension] + T[0][dimension];

      iterIls++;
    }

    if(better_cost < best_cost)
    {
      best_cost = better_cost;
      best_solution = better_solution;
    }
  }

  auto timerEnd = chrono::system_clock::now();
  //GILS-RVND
  constructReopt(W, T, C, best_solution); //atualiza estruturas auxiliares da best solution

  chrono::duration<double> timerGilsRVND = timerEnd - timerStart;

  std::cout << "TIME: " << timerGilsRVND.count() << "\n";
  printSolution(best_solution, best_cost, W, T, C);

  return 0;
}
