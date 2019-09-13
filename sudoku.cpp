/* This C++ program solves a Sudoku puzzle using Recursive Backtracking.
   Unsolved puzzle is accepted from user as an input file where unassigned (row,col) is marked as '0'. 
   It has 1 class - Sudoku, which contains 7 public member functions.The
   main() creates the object and calls the appropriate member functions to carry out the task*/


   //Written by: Ankur Mukherjee
   //Program: MSFE 2020, University of Illinois, Urbana-Champaign
   //Date: 13th September, 2019

#include<iostream> 
#include<fstream>
#include<string>

using namespace std;

const int UNASSIGNED = 0;
const int N = 9;

class Sudoku
{

public:

	int puzzle[9][9];
	// Public member function that (recursively) implements the brute-force 
	// search for possible solutions to the incomplete Sudoku puzzle
	bool Solve(int puzzle[N][N])
	{
		int row, col;

		
		if (!empty(puzzle, row, col))
			return true; 

			//looping down from 9 to 1 
		for (int num = 9; num >= 1; num--)
		{
			
			if (safeToAssign(puzzle, row, col, num))
			{
				
				puzzle[row][col] = num;

				
				if (Solve(puzzle))
					return true;

				// no solution, reassign to zero and start again 
				puzzle[row][col] = UNASSIGNED;
			}
		}
		return false; // backtracking 
	}

/* Traverses the puzzle to find an entry that is
still not assigned. If found, the reference
parameters row, col will be set the location
that is unassigned, and true is returned.
If no unassigned entries remain, false is returned. */
	bool empty(int puzzle[N][N],
		int& row, int& col)
	{
		for (row = 0; row < N; row++)
			for (col = 0; col < N; col++)
				if (puzzle[row][col] == UNASSIGNED)
					return true;
		return false;
	}

	//Public member function that checks if the named row is valid
	bool row_valid(int puzzle[N][N], int row, int num)
	{
		for (int col = 0; col < N; col++)
			if (puzzle[row][col] == num)
				return true;
		return false;
	}

	//Public member function that checks if the named col is valid
	bool col_valid(int puzzle[N][N], int col, int num)
	{
		for (int row = 0; row < N; row++)
			if (puzzle[row][col] == num)
				return true;
		return false;
	}

	//Public member function that checks if the named 3X3 grid is valid
	bool block_valid(int puzzle[N][N], int boxStartRow,
		int boxStartCol, int num)
	{
		for (int row = 0; row < 3; row++)
			for (int col = 0; col < 3; col++)
				if (puzzle[row + boxStartRow]
					[col + boxStartCol] == num)
					return true;
		return false;
	}

	//Test to check if 'num' can be assigned to row,col
	bool safeToAssign(int puzzle[N][N], int row,
		int col, int num)
	{
		
		return !row_valid(puzzle, row, num) &&
			!col_valid(puzzle, col, num) &&
			!block_valid(puzzle, row - row % 3,
				col - col % 3, num) &&
			puzzle[row][col] == UNASSIGNED;
	}

	//Public member function that prints the puzzle
	void print_puzzle(int puzzle[N][N])
	{
		for (int row = 0; row < N; row++)
		{
			for (int col = 0; col < N; col++)
				cout << puzzle[row][col] << " ";
			cout << endl;
		}
	}

	
};


void main()
{
	 
	string path;

	cout << "Entire path for the input file: " << endl;

	int puzzle[9][9];
	cin >> path;
	//path = C:\Users\ankur\source\repos\sudoku_1sttry\input1.txt

	ifstream infile(path);

	cout << "Unsolved Sudoku Puzzle: " << endl;
	//cout << infile.rdbuf();   //Using read buffer function


	if (infile.is_open())
	{
		//traversing till the end of the file
		while (!infile.eof())
		{
			for (int i = 0; i < 9; i++)
			{
				for (int j = 0; j < 9; j++)
				{
					//Storing the sudoku puzzle in a 2D array
					infile >> puzzle[i][j];
					
				}
				
			}
		}
	}

	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < 9; j++)
		{
			cout << puzzle[i][j] << " ";

		}
		cout << endl;
	}

	cout << "\nSolution:"<<endl;
	Sudoku x;

	
	if (x.Solve(puzzle) == true)
		x.print_puzzle(puzzle);
	else
		cout << "No solution exists";


}


