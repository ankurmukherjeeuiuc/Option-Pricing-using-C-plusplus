/* This C++ program prints one possible solution for N Queens problem using Recursive Backtracking.
   N is accepted from the user. It has 1 class - Board, which contains 4 public member functions.The
   main() creates the object and calls the appropriate member functions to carry out the task*/


   //Written by: Ankur Mukherjee
   //Program: MSFE 2020, University of Illinois, Urbana-Champaign
   //Date: 6th September, 2019


#include<iostream>
using namespace std;

class Board {


public:

	// member function: prints the board position
	void print_Board(int** board, int N) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++)
				cout << board[i][j] << " ";
			cout << endl;
		}
	}

	// member function:  returns 'false' if
	// the (row, col) position is not safe.

	bool is_this_position_safe(int** board, int row, int col, int N) {
		for (int i = 0; i < col; i++) //check whether there is queen in the left or not
			if (board[row][i])
				return false;
		for (int i = row, j = col; i >= 0 && j >= 0; i--, j--)
			if (board[i][j]) //check whether there is queen in the left upper diagonal or not
				return false;
		for (int i = row, j = col; j >= 0 && i < N; i++, j--)
			if (board[i][j]) //check whether there is queen in the left lower diagonal or not
				return false;
		return true;
	}

	//member function : recursive backtracking
	bool solve(int** board, int col, int N) {
		if (col >= N) //when N queens are placed successfully
			return true;
		for (int i = 0; i < N; i++) { /* Check if queen can be placed on
			board[i][col] */
			if (is_this_position_safe(board, i, col, N)) {
				board[i][col] = 1; /* Place this queen in board[i][col] */
				if (solve(board, col + 1, N)) //Go for the other columns recursively
					return true;
				board[i][col] = 0; //When no place is vacant remove that queen- Backtrack
			}
		}
		return false; //when no possible order is found
	}


	// Solves the n-Queens problem by (recursive) backtracking
	bool nQueens(int** board, int n) {

		if (solve(board, 0, n) == false) { //starting from 0th column
			cout << "Solution does not exist";
			return false;
		}
		print_Board(board, n);
		return true;
	}
};


void main(int argc, char* argv[]) {

	int N;
	cout << "Enter n:" << endl;
	cin >> N;

	//pointer to pointer for dynamic array initialization
	int** board = new int* [N];


	for (int i = 0; i < N; i++)
		board[i] = new int[N];

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			board[i][j] = 0; //set all elements to 0
		}
	}

	Board x;
	x.nQueens(board, N);
	

	}