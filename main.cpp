#include <stdio.h>
#include <iostream> 
#include <vector>
#include <cmath>

//////////////////////////////////////
//@brief triangle fuzzy number set structure
// a, b, c - left, middle and right values
struct TriangleFuzzySet
{
	double a;
	double b;
	double c;

	//TriangleFuzySet(double a0, double b0, double c0)
	//	:a(a0), b(b0), c(c0)
	//{};
	double MembershipVal(double x)
	{
		if(x <= a || x >= c)
			return 0;
		else if(x > a && x <= b) 
			return (x - a)/(b - a);
		else return -(x - b)/(c - b) + 1;

		return 0 ;		
	}

	double Norm()
	{
		const int N = 1000;
		double h = (c - a)/N;
		double s = 0;
		for(double x = a; x <= c; x += h)
		{
			s += MembershipVal(x)*h;
		}

		return s;
	}
};

//////////////////////////////////////
//@brief
void LoadInputSets(std::vector<TriangleFuzzySet> &A, const char *fname = "in.txt")
{
	FILE *f;
	int N;
	f = fopen(fname, "r");
	fscanf(f, "%d", &N);
	A.resize(N);

	for(int i = 0; i < N; ++i)
		fscanf(f, "%lf%lf%lf", &A[i].a, &A[i].b, &A[i].c);

	fclose(f);
}

//////////////////////////////////////
//@brief Hemming distance between Triangle Fuzzy Sets A and B
double Hemming(TriangleFuzzySet &A, TriangleFuzzySet &B)
{
	const int N = 1000;

	double x_min = std::min(A.a, B.a);
	double x_max = std::max(A.c, B.c);
	double h = (x_max - x_min)/(double)N;

	double s = 0;
	for(double x = x_min; x <= x_max; x += h)
	{
		//double t = abs(A.MembershipVal(x) - B.MembershipVal(x));
		s += abs(A.MembershipVal(x) - B.MembershipVal(x))*h;
	}
	double l = A.Norm() + B.Norm();

	return s/l;
}

//////////////////////////////////////
//@brief
std::vector<std::vector<double> > HemmingMatrix(std::vector<TriangleFuzzySet> &A)
{
	int size = (int) A.size();
	std::vector<std::vector<double> > R(size);
	for(int i = 0; i < size; ++i)
	{
		R[i].resize(size);
		for(int j = 0; j < size; ++j)
		{
			if(i < j)
			{
				R[i][j] = Hemming(A[i], A[j]);
			}
			else
			{
				if(i == j)
				{
					R[i][j] = 0;
				}
				else
				{
					R[i][j] = R[j][i];
				}
			}
		}
	}
	return R;
}


std::vector<std::vector<double> > LoadInputMatrix(/*std::vector<std::vector<double> > &in*/)
{
	std::vector<std::vector<double> > in;
	char *fname = "in";
	//Чтение из файла исходной информации
	int size = 0;
	FILE *f = fopen(fname, "r+");
	fscanf(f, "%d", &size);
	in.resize(size);
	for(int i = 0; i < size; ++i)
	{
		in[i].resize(size);
		for(int j = 0; j < size; ++j)
		{
			fscanf(f, "%lf", &in[i][j]);
		}
	}

	fclose(f);
	return in;
}

void PrintMatrix(std::vector<std::vector<double> > &in, const char *str){
	printf(str);
	for(int i = 0; i < (int) in.size(); ++i)
	{
		for(int j = 0; j < (int) in.size(); ++j)
		{
			printf("%1.3lf\t", in[i][j]);
		}
		printf("\n");
	}
}


//////////////////////////////////////
//@brief
template <class T> void PrintMatrixT(std::vector<std::vector<T> > &M, const char *str)
{
	std::cout << str;
	for(int i = 0; i < M.size(); ++i)
	{
		for(int j = 0; j < M.size(); ++j)
		{
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

double Hemming(std::vector<double> &A, std::vector<double> &B)
{
	double l = 0;
	for(int i = 0; i < (int) A.size(); ++i)
	{
		l += fabs(A[i] - B[i]);
	}
	return l/A.size();
}

std::vector<std::vector<double> > HemmingMatrix(std::vector<std::vector<double> > &in)
{
	int size = (int) in.size();
	std::vector<std::vector<double> > R(size);
	for(int i = 0; i < size; ++i)
	{
		R[i].resize(size);
		for(int j = 0; j < size; ++j)
		{
			if(i < j)
			{
				R[i][j] = Hemming(in[i], in[j]);
			}
			else
			{
				if(i == j)
				{
					R[i][j] = 0;
				}
				else
				{
					R[i][j] = R[j][i];
				}
			}
		}
	}
	return R;
}

//////////////////////////////////////
//@brief
std::vector<std::vector<double> > MinMaxComposition(const std::vector<std::vector<double> > &R0, const std::vector<std::vector<double> > &R)
{
	std::vector<std::vector<double> > R_2(R.size());
	int N = R.size();
	for(int i = 0; i < N; ++i)
	{
		R_2[i].resize(N);
		for(int j = 0; j < N; ++j)
		{
			double Rij_min = 1;
			for(int k = 0; k < N; ++k)
			{
				double m = std::max(R0[i][k], R[k][j]);
				if(m < Rij_min) 
					Rij_min = m;
			}
			R_2[i][j] = Rij_min;
		}			
	}

	return R_2;
}

//////////////////////////////////////
//@brief compare of two double valued symetric matrix
bool VecCmp(const std::vector<std::vector<double> > &R1, const std::vector<std::vector<double> > &R2)
{
	const double eps = 1e-3;
	bool flag = true;

	for(int i = 0; i < R1.size(); ++i){
		for(int j = 0; j < R1.size(); ++j)
			if(std::abs(R1[i][j]-R2[i][j]) > eps) 
			{
				flag = false;
				printf ("Vectors are different in: i=%d, j=%d\n", i, j);
				break;
			}
			if(flag == false) break;
	}
			return flag;
}

////////////////////////////////////////
//@ MinMax transitive closure computing
std::vector<std::vector<double> > MinMaxTransitiveClosure(const std::vector<std::vector<double> > &R)
{
	int N = R.size();
	std::vector<std::vector<double> > R_TrCl_1(N), R_TrCl(N);
	int i = 1;
	R_TrCl = R;
	do
	{
		R_TrCl_1 = R_TrCl;
		R_TrCl = MinMaxComposition(R_TrCl, R);
		char title[255];
		sprintf(title, "\nR^%d\n", ++i);
		PrintMatrix(R_TrCl, title);

	}
	while(!VecCmp(R_TrCl, R_TrCl_1));

	return R_TrCl;
}

//////////////////////////////////////
//@brief
void RelationNegation(std::vector<std::vector<double> > &R)
{
	for(int i = 0; i < R.size(); ++i)
		for(int j = 0; j < R.size(); ++j)
			R[i][j] = 1 - R[i][j];
}

//////////////////////////////////////
//@brief
std::vector<std::vector<int> > AlphaLevel(const std::vector<std::vector<double> > &R, double alpha)
{
	int N = R.size();	
	std::vector<std::vector<int> >  Falpha(N);

	for(int i = 0; i < N; ++i)
	{
		Falpha[i].resize(N, 0);
		for(int j = 0; j < N; ++j)
			if(R[i][j] >= alpha) Falpha[i][j] = 1;

	}

	return Falpha;
}

//////////////////////////////////////
//@brief
void RowColSwap(std::vector<std::vector<int> > &R, int i, int j)
{
	std::swap(R[i], R[j]);
	for(int k = 0; k < R.size(); ++k)
		std::swap(R[k][i], R[k][j]);
}

//////////////////////////////////////
//@brief function making block matrix from matrix R
void BlockMatrix(std::vector<std::vector<int> > &R, std::vector<int> &tr)
{
	int N = R.size();
	tr.resize(N);
	for(int i = 0; i < N; ++i)
		tr[i] = i;

	int i0 = 0, i1, row = 0;
	while(row < N - 2)
	{
		while((R[row][++i0] != 0) && (i0 < N-2));

		if(R[row][i0] == 0)
		{
			i1 = i0;
			while((i1 < N-1) && (R[row][++i1] != 1));

			if(R[row][i1] == 1)
			{
				RowColSwap(R, i0, i1);
				std::swap(tr[i0], tr[i1]);
			}
			else row = i0;

		} 
		else row = i0;
	}
}

//////////////////////////////////////
//@brief
void PrintClasses(std::vector<std::vector<int> > &R, std::vector<int> &tr)
{
	int N = R.size();
	int i = 0, ind = 0, k = 0;
	//printf("\nPermutations of columns:\n");
	for(int i = 0; i < tr.size(); ++i)
		printf("%d ", tr[i]);
	printf("\n");
	printf("\n\tClasses:\n");
	while(i < N)
	{
		printf("L%d = { ", ind++);
		//k = i;
		while(k < N && R[i][k])
		{
			printf("%d ", tr[k]);
			++k;
		}
		printf("}\n");
		i = k;
	}
}

//////////////////////////////////////
//@brief                                                                                                                                                                                         100
double GetNextMin(std::vector<std::vector<double> > &R, double cur_min)
{
	const double eps = 1e-9;
	double min = 1;
	if(cur_min == 1) return 2;

	for(int i = 0; i < R.size(); ++i)
		for(int j = 0; j <= i; ++j)
			if((R[i][j] > cur_min + eps) && (R[i][j] < min))
				min = R[i][j];
	return min;
}

//////////////////////////////////////
//@brief
void FuzzyClassification(std::vector<std::vector<double> > &R)
{
	std::vector<int> tr;
	std::vector<std::vector<int> > Ra;
	double alpha = 0;
	for(alpha = GetNextMin(R, alpha); alpha <= 1; alpha  = GetNextMin(R, alpha))
	{
		Ra = AlphaLevel(R, alpha);
		printf("\n\talpha = %lf", alpha);
		//PrintMatrixT<int>(Ra, "\nR:\n");
		BlockMatrix(Ra, tr);
		PrintMatrixT<int>(Ra, "\nBlock R :\n");
		PrintClasses(Ra, tr);
	}
}


/////////////////////////////////
//@brief main
int main (void)
{
	std::vector<std::vector<double> > in;
	std::vector<TriangleFuzzySet> A;
	std::vector<std::vector<double> > R;
	std::vector<std::vector<double> > R_cl;
	int i;

	do
	{

		printf("Choose:\n1. 8 elemens finite fuzzy set.\n2. Triangular fuzzy sets.\n");
		scanf("%d", &i);
		switch(i)
		{
		case 1:
			in = LoadInputMatrix();
			R = HemmingMatrix(in);
			PrintMatrix(in, "\nInput Matrix\n");
			break;
		case 2:
			LoadInputSets(A);
			R = HemmingMatrix(A);
			break;
		default: 
			printf("Wrong value\n");
		}
	}while(i != 1 && i !=2);

	PrintMatrix(R, "\nMatrix R:\n");

	R_cl = MinMaxTransitiveClosure(R);
	PrintMatrix(R_cl, "\nTransitive closure of R:\n");

	RelationNegation(R_cl);
	PrintMatrix(R_cl, "\nEquality relation :\n");

	FuzzyClassification(R_cl);

	getchar();

	getchar();
	return 0;
}
