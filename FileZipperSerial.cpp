#include <iostream>
#include <fstream>
#include <omp.h>
using namespace std;
#pragma pack(2)
class BMPHeader
{
	private:
		char HeaderField[2];
		unsigned int size;
		unsigned int garbage;
		unsigned int imageOffset;
	public:
		unsigned int getImageOffset(){return imageOffset;}
};
class DIBHeader
{
	private:
		unsigned int hsize;
		int width;
		int height;
		unsigned short int colorplanes;
		unsigned short int bpp;
		unsigned int compression;
		unsigned int size;
		unsigned int hor_res;
		unsigned int ver_res;
		unsigned int colorDepth;
		unsigned int impColor;
	public:
		unsigned int getWidth(){return width;}
};
class CharInfo
{
	private:
		int freq;
		char c;
		unsigned short int BinaryCode;
		int BinaryLength;
	public:
		//Consructors
		CharInfo(){freq = 0; c = '\0'; BinaryCode = 0b00000000; BinaryLength = 0;}
		CharInfo(int f, char c){this->freq = f; this->c = c; BinaryCode = 0b00000000; BinaryLength = 0;}
		CharInfo(int f, char c, unsigned short int bin, int len){this->freq = f; this->c = c; BinaryCode = bin; BinaryLength = len;}
		//Setters
		void setBinaryCode(unsigned short int val){BinaryCode = val;}
		void setBinaryLength(int val){BinaryLength = val;}
		void setChar(char c){this->c = c;}
		void setfreq(int f){this->freq = f;}
		//Getters
		bool getBinaryCode(int index){return (BinaryCode & (1<<index))>>index;}
		int getBinaryLength(){return BinaryLength;}
		char getChar(){return c;}
		int getFreq(){return freq;}
		//Overloaded Operators
		bool operator > (CharInfo &x) {return freq > x.freq;}
		bool operator < (CharInfo &x) {return freq < x.freq;}
		void operator ++ () {++freq;}
		void operator = (CharInfo &x) {this->freq = x.freq; this->c = x.c; this->BinaryLength = x.BinaryLength; BinaryCode = x.BinaryCode;}
		CharInfo operator + (CharInfo &x) {return CharInfo(this->freq + x.freq,'§');}
		bool operator == (char c) {return this->c == c;}
		//Member Functions
};
template<typename T>
class Node
{
	public:
		T data;
		Node<T> *left, *right;
	public:
		//Constructors
		Node(){left = right = NULL;}
		Node(T data){left = right = NULL; this->data = data;}
		//Destructors
		~Node(){left = right = NULL;}
		//Overloaded Operators
		bool operator > (Node<T> &x) {return data > x.data;}
		bool operator < (Node<T> &x) {return data < x.data;}
		void operator = (Node<T> *x) {*data = x->data;}
		//Member Functions
};
template<typename T>
class LinkedList
{
	Node<T> *head;
	public:
		//Constructors
		LinkedList():head(NULL){}
		LinkedList(T d){head = new Node<T>(d);}
		//Destructor
		~LinkedList(){DeleteList();}
		//Member Functions
		void DeleteList()
		{
			if (head == NULL) return;
			else
			{
				while(head != NULL)
				{
					Node<T>* temp = head;
					head = head->right;
					delete temp;
				}
			}
		}
		void Insert(T d)
		{
			Node<T>* newNode = new Node<T>(d);
			newNode->right = head;
			head = newNode;
		}
		Node<T>* ReturnHead(){return head;}
		void InsertChar(char ch, int value = 1)
		{
			Node<T>* temp = head;
			while(temp != NULL)
			{
				if (temp->data == ch)
				{
					temp->data.setfreq(temp->data.getFreq()+value);
					return;
				}
				temp = temp->right;
			}
			Insert(CharInfo(1,ch));
			return;
		}
		Node<T>* Pop()
		{
			Node<T>* temp = head;
			head = head->right;
			temp->right = NULL;
			return temp;
		}
		int SizeofList()
		{
			int x = 0;
			Node<T> *temp = head;
			while(temp != NULL)
			{
				temp = temp->right;
				x++;
			}
			return x;
		}
};
template<typename T>
class MinHeap
{
	Node<T> **arr;
	int capacity;
	int top;
	int parent(int i) {return (i-1)/2;}
	int left(int i) {return 2*i + 1;}
	int right(int i) {return 2*i + 2;}
	public:
		//Constructors
		MinHeap(){arr = NULL;}
		MinHeap(int n)
		{
			capacity = n;
			top = 0;
			arr = new Node<T>*[capacity];
		}
		//Destructor
		~MinHeap()
		{
			capacity = 0;
			top = 0;
			delete [] arr;
		}
		//Member Functions
		int Capacity(){return capacity;}
		void Insert(Node<T> *d)
		{
			if (top == capacity) return;
			int i = top++;
			arr[i] = d;
			while(i != 0 && *arr[parent(i)] > *arr[i])
			{
				Node<T> *x = arr[parent(i)];
				arr[parent(i)] = arr[i];
				arr[i] = x;
				i = parent(i);
			}
		}
		void InsertList(LinkedList<T>& List) {for(int i = 0 ; i < capacity ; i++) Insert(List.Pop());}
		Node<T>* ExtractRoot()
		{
			if(top > 0)
			{
				if(top == 1)
				{
					top--;
					return arr[0];
				}
				Node<T> *root = arr[0];
				arr[0] = arr[top - 1];
				top--;
				MinHeapify(0);
				return root;
			}
		}
		void MinHeapify(int i)
		{
			int leftIndex = left(i);
			int rightIndex = right(i);
			int smallest = i;
			if(leftIndex < top && *arr[leftIndex] < *arr[smallest]) smallest = leftIndex;
			if(rightIndex < top && *arr[rightIndex] < *arr[smallest]) smallest = rightIndex;
			if(smallest != i) 
			{
				Node<T> *x = arr[i];
				arr[i] = arr[smallest];
				arr[smallest] = x;
				MinHeapify(smallest);
 			}
		}
		Node<T>* CreateHuffmanTree()
		{
			Node<T> *left, *right, *root;
			while(top != 1)
			{
				left = ExtractRoot();
				right = ExtractRoot();
				root = new Node<T>(left->data+right->data);
				root->left = left;
				root->right = right;
				Insert(root);
			}
			return ExtractRoot();
		}
};
template<typename T>
class HuffmanTree
{
	Node<T>* root;
	public:
		//Constructors
		HuffmanTree(){root = NULL;}
		HuffmanTree(Node<T>* n){root = n;}
		//Destructor
		~HuffmanTree(){DeleteHuffmanTree();}
		//Member Functions
		DeleteHuffmanTree(){DeleteHuffmanTree(root);root=NULL;}
		void DeleteHuffmanTree(Node<T>* r)
		{
			if(r == NULL)return;
			DeleteHuffmanTree(r->left);
			DeleteHuffmanTree(r->right);
			delete r;
		}
		//Encoding & Compression Functions
		void EncodeHuffmanCodes(fstream& file2 , LinkedList<T>& L1, bool NOWRITE = 0, int SKIP = 0)
		{
			unsigned short int BinaryCode;
			EncodeHuffmanCodes(root,BinaryCode,0,file2,L1,NOWRITE);
			if(!NOWRITE)
			{
				file2 << root->data.getFreq();
				file2 << '§';
			}
		}
		void EncodeHuffmanCodes(Node<T> *r , unsigned short int &t , int index , fstream& File2 , LinkedList<T>& L, bool NOWRITE = 0)
		{
			if(r->left)
			{
				t &= ~(1<<index);
				EncodeHuffmanCodes(r->left, t, index + 1, File2 , L, NOWRITE);
			}
			if(r->right)
			{
				t |= (1<<index);
				EncodeHuffmanCodes(r->right, t, index + 1, File2 , L, NOWRITE);
			}
			if(!r->left && !r->right)
			{
				r->data.setBinaryLength(index);
				r->data.setBinaryCode(t);
				if(!NOWRITE)
				{
					File2 << index;
					File2 << '§';
					File2 << r->data.getChar();
					File2 << t;
					File2 << '§';
				}
				L.Insert(r->data);
			}
		}
		void CompressFile(fstream& File1 , fstream& File2 , LinkedList<T>& L, int Width = 0, int SKIP = 0)
		{
			unsigned char x = 0;
			char c;
			int count = 0;
			int POSCounter = 0;
			while(File1.get(c))
			{
				POSCounter++;
				Node<T>* temp = L.ReturnHead();
				while(temp->data.getChar() != c && temp->right) temp = temp->right;
				if(temp->data.getChar() == c)
				{
					for(int i = 0 ; i < temp->data.getBinaryLength() ; i++)
					{	
						if(count < 7)
						{
							if(temp->data.getBinaryCode(i)) x++;
							x <<= 1;
							count++;
						}
						else if (count == 7)
						{
							if(temp->data.getBinaryCode(i)) x++;
							count = 0;
							File2 << x;
							x = 0;
						}
					}
				}
				if(POSCounter == Width*3)
				{
					POSCounter = 0;
					File1.seekg(SKIP, ios::cur);
				}
			}
			for (int i = 0 ; i < 7 - count ; i++) x <<= 1;
			File2 << x;
		}
		//Decoding and Decompression Functions
		void DecodeHuffmanCodes(fstream& file2)
		{
			int n;
			char EOC;
			file2 >> n;
			file2 >> EOC;
			DeleteHuffmanTree();
			root = DecodeHuffmanCodes(file2,root,n);
		}
		Node<T>* DecodeHuffmanCodes(fstream& file2, Node<T> *r, int size)
		{
			r = new Node<T>(CharInfo(0,'§'));
			for (int i = 0 ; i < size ; i++)
			{
				char c, EOC;
				int length;
				unsigned short int code;
				file2 >> length;
				file2 >> EOC;
				file2.get(c);
				file2 >> code;
				file2 >> EOC;
				Node<T>* temp = r;
				for(int j = 0 ; j < length ; j++)
				{
					if(((code & (1<<j))>>j))
					{
						if(!temp->right) temp->right = new Node<T>(CharInfo(0,'§',code,length));
						temp = temp->right;
					}
					else
					{
						if(!temp->left) temp->left = new Node<T>(CharInfo(0,'§',code,length));
						temp = temp->left;
					}
				}
				temp->data.setChar(c);
			}
			return r;
		}
		void DecompressFile(fstream& file2 , fstream& file3, int Width = 0, int SKIP = 0)
		{
			unsigned char x;
			Node<T>* r = root;
			char EOC;
			int i = 7;
			long long int k, POSCounter = 0;
			file2 >> k;
			file2 >> EOC;
			while(true)
			{
				x = file2.get();
				while(true)
				{
					if(k <= 0) return;
					if(!r->left && !r->right)
					{
						POSCounter++;
						file3 << r->data.getChar();
						k--;
						r = root;
						if(POSCounter == Width*3)
						{
							POSCounter = 0;
							EOC = r->data.getChar();
							for(int z = 0 ; z < SKIP ; z++) file3.write((char*)&EOC, 1);
						}
						continue;
					}
					if(i == -1)
					{
						i = 7;
						break;
					}
					if((x & (1<<i))>>i) r = r->right;
					else r = r->left;
					i--;
				}
			}
		}
};
void ZipFile()//Main Zip Function
{
	LinkedList<CharInfo> List;
	fstream ReadFile, WriteFile;
	string FileName = "";
	char c;

	cout << "Enter Name of the File to Zip: ";
	cin >> FileName;
	
	//Opening File To ZIP
	ReadFile.open(FileName.c_str() , ios::in | ios::binary);
	if(!ReadFile.is_open()) {cout << "No Such File Exists!";return;}
	
	if(!FileName.compare(FileName.find('.'),4,".txt"))
	{
		//Creating Compressed File
		for(int i = 0 ; i < 4 ; i++) FileName.pop_back();
		FileName += "Cmp.txt";
		WriteFile.open(FileName.c_str(), ios::out | ios::binary);
		
		//Reading All Characters from File
		while(ReadFile.get(c)) List.InsertChar(c);
		cout << "Compressing File...\n";
		
		//Putting all unique characters into Priority Queue
		MinHeap<CharInfo> PriorityQueue(List.SizeofList());
		PriorityQueue.InsertList(List);
		WriteFile << PriorityQueue.Capacity() << "§";
	
		//Creating Huffman Tree then Encoding and Compressing Data
		HuffmanTree<CharInfo> Tree(PriorityQueue.CreateHuffmanTree());
		Tree.EncodeHuffmanCodes(WriteFile, List);
		ReadFile.clear();
		ReadFile.seekg(0);
		Tree.CompressFile(ReadFile, WriteFile, List);
		WriteFile.close();
		cout << "File Successfully Zipped with Name: " << FileName << "!\n";
	}
	else if(!FileName.compare(FileName.find('.'),4,".bmp"))
	{
		//Creating Compressed File
		for(int i = 0 ; i < 4 ; i++) FileName.pop_back();
		FileName += "Cmp.bmp";
		WriteFile.open(FileName.c_str(), ios::out | ios::binary);
		
		//Reading Image File Headers
		BMPHeader bheader;
		DIBHeader dheader;
		ReadFile.read((char*)&bheader, sizeof(bheader));
		if(bheader.getImageOffset() != 54)
		{
			cout << "Unsupported Format of BMP! 24-bit is supported!" << endl;
			return;
		}
		ReadFile.read((char*)&dheader,sizeof(dheader));
		cout << "Compressing File...\n";
		//Reading All Characters from File
		double t1 = omp_get_wtime();
		ReadFile.seekg(bheader.getImageOffset(), ios::beg);
		int count = 0;
		int Padding = ((4 - (dheader.getWidth() * 3)%4)%4);
		while(ReadFile.get(c))
		{
			List.InsertChar(c);
			count++;
			if(count == dheader.getWidth()*3)
			{
				count = 0;
				ReadFile.seekg(Padding, ios::cur);
			}
		}
		
		//Putting all unique characters into Priority Queue
		MinHeap<CharInfo> PriorityQueue(List.SizeofList());
		PriorityQueue.InsertList(List);
		WriteFile.write((char*)&bheader,sizeof(bheader));
		WriteFile.write((char*)&dheader,sizeof(dheader));
	
		//Creating Huffman Tree then Encoding and Compressing Data
		WriteFile.seekp(bheader.getImageOffset(), ios::beg);
		WriteFile << PriorityQueue.Capacity() << "§";
		HuffmanTree<CharInfo> Tree(PriorityQueue.CreateHuffmanTree());
		Tree.EncodeHuffmanCodes(WriteFile, List);
		ReadFile.clear();
		ReadFile.seekg(bheader.getImageOffset(), ios::beg);
		Tree.CompressFile(ReadFile, WriteFile, List, dheader.getWidth(), Padding);
		double t2 = omp_get_wtime();
	
		WriteFile.close();
		cout << "File Successfully Zipped with Name: " << FileName << "!\n";
		cout << "Time taken: " << t2-t1 << "s\n";
	}
	else cout << "Exception Thrown! Incorrect or Unsupported file type entered!";
	ReadFile.close();
}
void UnzipFile()//Main Unzip Function
{
	fstream ReadFile, WriteFile;
	string FileName = "";
	
	cout << "Enter Name of the File to Unzip: ";
	cin >> FileName;
	
	//Opening Compressed File
	ReadFile.open(FileName.c_str() , ios::in | ios::binary);
	if(!ReadFile.is_open()) {cout << "No Such File Exists!";return;}

	if(!FileName.compare(FileName.find('.'),4,".txt"))
	{
		//Creating Decompressed File
		for(int i = 0 ; i < 7 ; i++) FileName.pop_back();
		FileName += "Decmp.txt";
		WriteFile.open(FileName.c_str() , ios::out | ios::binary);
		cout << "Decompressing File...\n";
	
		//Recreating Huffman Tree and Decompressing File
		HuffmanTree<CharInfo> Tree;
		Tree.DecodeHuffmanCodes(ReadFile);
		Tree.DecompressFile(ReadFile, WriteFile);
	
		WriteFile.close();
	}
	else if(!FileName.compare(FileName.find('.'),4,".bmp"))
	{
		for(int i = 0 ; i < 7 ; i++) FileName.pop_back();
		FileName += "Decmp.bmp";
		WriteFile.open(FileName.c_str(), ios::out | ios::binary);
		
		BMPHeader bheader;
		DIBHeader dheader;
		cout << "Decompressing File...\n";
		
		ReadFile.read((char*)&bheader, sizeof(bheader));
		if(bheader.getImageOffset() != 54)
		{
			cout << "Unsupported Format of BMP! 24-bit is supported!" << endl;
			return;
		}
		ReadFile.read((char*)&dheader, sizeof(dheader));
		
		WriteFile.write((char*)&bheader, sizeof(bheader));
		WriteFile.write((char*)&dheader,sizeof(dheader));
		ReadFile.seekg(bheader.getImageOffset(), ios::beg);
		HuffmanTree<CharInfo> Tree;
		Tree.DecodeHuffmanCodes(ReadFile);
		int Padding = ((4 - (dheader.getWidth() * 3)%4)%4);
		Tree.DecompressFile(ReadFile, WriteFile, dheader.getWidth(), Padding);
		WriteFile.close();
	}
	ReadFile.close();
	cout << "File Successfully Unzipped with Name: " << FileName << "!\n";
}
void PrintFile()//Print Contents of File
{
	fstream ReadFile;
	string temp;
	
	cout << "Enter Name of the File to print: ";
	cin >> temp;
	
	//Opening File To Print
	temp += ".txt";
	const char *ReadFileTxt = temp.c_str();
	ReadFile.open(temp.c_str() , ios::in | ios::binary);
	if(!ReadFile.is_open()) {cout << "No Such File Exists!";return;}
	
	//Reading Contents of File
	cout << '\n';
	while(true)
	{
		char c;
		ReadFile.get(c);
		if(ReadFile.eof()) break;
		cout << c;
	}
	ReadFile.close();
}
int main()
{
	int choice;
	while(true)
	{	
		cout << "\t\tFAST NUCES TXT File Zipper\n" << "1.Zip File\n2.Unzip File\n3.Print File\n4.Exit\nChoice: ";
		cin >> choice;
		if(choice == 1) ZipFile();
		else if (choice == 2) UnzipFile();
		else if (choice == 3) PrintFile();
		else if (choice == 4) break;
		cin.get();
		cin.get();
		system("CLS");
	}	
	cout << "\nMade By Mohammad Yehya Hayati (K213309), Mahad Munir (K213388), Daniyal Naqvi (K213433)";
}
