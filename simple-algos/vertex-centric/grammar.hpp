#ifndef GRAMMAR_HPP
#define GRAMMAR_HPP

#include "globals.hpp"


class Grammar{
public:
	string grammarFilePath;
	// level-1: rule-ID, level-2: LHS [0]
	vector<vector<uint>> grammar1;
	// level-1: rule-ID, level-2: LHS [0], RHS [1]
	vector<vector<uint>> grammar2;
	// SF:: Let's say A --> BC is the grammar rule. 
	// Then first level of the 2D list will be the number when it's found among it's kind. If we have seen this kind of
	// rule before once then it's index will be 2. Second index will be the grammar rules.
	// For example grammar3[3][0] = A, grammar3[3][1] = B, grammar3[3][2] = C
	// level-1: rule-ID, level-2: LHS [0], RHS-1 [1], RHS-2 [2]
	vector<vector<uint>> grammar3;
	// index: label (LHS), value: in the grammar or not
	vector<bool> grammar1index;
	// consider A --> B, index: B, value: all possible As
	vector<vector<uint>> grammar2index;
	// consider A --> BC, index: BC (B*num_symbols+C), value: all possible As
	vector<vector<uint>> grammar3index;
	// consider A --> BC, index: B, value: all possible pairs of (A, C)
	vector<vector<pair<uint, uint>>> grammar3indexLeft;
	// consider A --> BC, index: C, value: all possible pairs of (A, B)
	vector<vector<pair<uint, uint>>> grammar3indexRight;
    // rule index

	int labelSize;
	unordered_map<string, uint> hashSym; // 0 .. labelSize-1
	unordered_map<uint, string> hashSymRev;
	int leftLabelSize;
	vector<uint> leftLabels;
	unordered_set<uint> allLabels;
	

	// SF: Reads the grammar file

	Grammar(string grammarFilePath)
	{
		this->grammarFilePath = grammarFilePath;
		ifstream grammarFile;
		grammarFile.open(grammarFilePath);
		
		string line;
		string symbol;
		int numSym;
		vector<uint> symbols;

		stringstream ss;
		labelSize = 0;
		while(getline(grammarFile, line))
		{
			ss.clear();
			symbols.clear();
			ss << line;
			numSym = 0;
			while(ss >> symbol)
			{
				if(hashSym.count(symbol) <= 0)
				{
					hashSym[symbol] = labelSize;
					hashSymRev[labelSize] = symbol;
					labelSize++;
				}
				numSym++;
				symbols.push_back(hashSym[symbol]); // pushing the index (here labelSize) of the symbol instead of the symbol itself
			}
			if(numSym == 1)
				// grammar1[0][0] will return the index of the symbol. Then we can use it to get the 
				// actual symbol from the hashSymRev
				grammar1.push_back(symbols); 
			else if (numSym == 2)
				grammar2.push_back(symbols);
			else if(numSym == 3)
				grammar3.push_back(symbols);
			else
			{
				cout << "An error happened during parsing the grammar!\n";
				exit(0);
			}
		}
		
		grammar1index.resize(labelSize);
		grammar2index.resize(labelSize);
		grammar3index.resize(labelSize*labelSize);
		grammar3indexLeft.resize(labelSize);
		grammar3indexRight.resize(labelSize);
		
		for(int i=0; i<labelSize; i++)
			grammar1index[i] = false;
		for(int i=0; i<grammar1.size(); i++) {
			grammar1index[grammar1[i][0]] = true;
			allLabels.insert(grammar1[i][0]);
		}
			
		for(int i=0; i<grammar2.size(); i++) {
			grammar2index[grammar2[i][1]].push_back(grammar2[i][0]);
         
			allLabels.insert(grammar2[i][0]);
			allLabels.insert(grammar2[i][1]);
		}
		for(int i=0; i<grammar3.size(); i++) {
			grammar3index[grammar3[i][1]*labelSize+grammar3[i][2]].push_back(grammar3[i][0]);

			allLabels.insert(grammar3[i][0]);
			allLabels.insert(grammar3[i][1]);
			allLabels.insert(grammar3[i][2]);
		}
		
		for(int i=0; i<grammar3.size(); i++)
			grammar3indexLeft[grammar3[i][1]].push_back(make_pair(grammar3[i][2], grammar3[i][0]));
		for(int i=0; i<grammar3.size(); i++)
			grammar3indexRight[grammar3[i][2]].push_back(make_pair(grammar3[i][1], grammar3[i][0]));
			
		
		unordered_set<uint> leftSym;
		for(int i=0; i<grammar2.size(); i++)
			leftSym.insert(grammar2[i][0]);
		for(int i=0; i<grammar3.size(); i++)
			leftSym.insert(grammar3[i][0]);
		
		for (auto it = leftSym.begin(); it != leftSym.end(); ++it) {
			leftLabels.push_back( (*it) );
		}
		leftLabelSize = leftLabels.size();
		
		cout << "Left Labels:\n";
		for(int i=0; i<leftLabelSize; i++)
			cout << hashSymRev[leftLabels[i]] << " ";
		cout << endl;
		
		
			 
		
		grammarFile.close();
	}
	
	
	
	inline bool rule1(uint symbol)
	{
		return grammar1index[symbol];
	}
	
	inline vector<uint> rule2(uint symbol)
	{
		return grammar2index[symbol];
	}
	
	inline vector<uint> rule3(uint symbol1, uint symbol2)
	{
		return grammar3index[symbol1*this->labelSize + symbol2];
	}
	
	inline vector<pair<uint, uint>> rule3left(uint symbol)
	{
		return grammar3indexLeft[symbol];
	}
	
	inline vector<pair<uint, uint>> rule3right(uint symbol)
	{
		return grammar3indexRight[symbol];
	}

    void printGrammarIndex3()
    {
        cout << "print1" <<endl;
        for (uint i = 0; i < this->grammar3.size(); i++)
        {
            cout << "print2" << endl;
            vector<uint> leftLabels = this->rule3(this->grammar3[i][1], this->grammar3[i][2]);
            for (uint j = 0; j < leftLabels.size(); j++)
            {
                cout << "print3" << endl;
                cout <<"(" << this->grammar3[i][1] << "," << this->grammar3[i][2] << ") ->" << leftLabels[j] << ",  ";
            }
            cout << endl;
        }
    }	


    void printGrammarIndex3_2()
    {
        cout << "print1" << this->labelSize <<endl;
        for (uint i = 0; i < this->grammar3index.size(); i++)
        {
            cout << "print2" << endl;
            vector<uint> leftLabels = this->grammar3index[i];
            for (uint j = 0; j < leftLabels.size(); j++)
            {
                cout << "print3" << endl;
                cout <<"(" << i << ") ->" << leftLabels[j] << ",  ";
            }
            cout << endl;
        }
    }

};


#endif
