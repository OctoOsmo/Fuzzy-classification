#include <stdio.h>
#include <vector>
#include <cmath>

std::vector<std::vector<double> > LoadInputMatrix(/*std::vector<std::vector<double> > &in*/)
{
    std::vector<std::vector<double> > in;
    //Чтение из файла исходной информации
    int size = 0;
    FILE *f = fopen("in.txt", "r+");
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
    return in;
}

void PrintMatrix(std::vector<std::vector<double> > &in, const char *str){
    printf(str);
    for(int i = 0; i < (int) in.size(); ++i)
    {
        for(int j = 0; j < (int) in.size(); ++j)
        {
            printf("%lf ", in[i][j]);
        }
        printf("\n");
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

bool VecCmp(const std::vector<std::vector<double> > &R1, const std::vector<std::vector<double> > &R2)
{
	bool flag = true;

	for(int i = 0; i < R1.size(); ++i)
		for(int j = 0; j < R1.size(); ++j)
			if(R1[i][j] != R2[i][j]) 
			{
				flag = false;
				break;
			}
	return flag;
}

std::vector<std::vector<double> > MinMaxTransitiveClosure(const std::vector<std::vector<double> > &R)
{
	int N = R.size();
	std::vector<std::vector<double> > R_TrCl_1(N), R_TrCl(N);

	R_TrCl = R;
	do
	{
		R_TrCl_1 = R_TrCl;
		R_TrCl = MinMaxComposition(R_TrCl, R);
	}
	while(VecCmp(R_TrCl, R_TrCl_1));
	
	return R_TrCl;
}

void RelationNegation(std::vector<std::vector<double> > &R)
{
	for(int i = 0; i < R.size(); ++i)
		for(int j = 0; j < R.size(); ++j)
			R[i][j] = 1 - R[i][j];
}

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

int GetNextOne(std::vector<int> v, int &k)
{
	int i = k;
	while(v[i++] == 0 && i <= v.size());

	if(i > v.size()) return -1;
	else return i;

}

void Swap(std::vector<std::vector<int> > &R, &k`
{
	int N = R.size();
	for(int i = 0; i < N; ++i)
		for(int j =0; j < N; ++j)

}

void DecompositionTree(std::vector<std::vector<int> > &R)
{
	int N = R.size();
	std::vector<int> transp(N);

	for(int i = 0; i < N; ++i)
		transp[i] = i;

	for(int i = 0; i < N; ++i)
	{
		for(int j = 0; j < N - i; ++j)
		{
			if(GetNextOne(R[i], j)

				
		}
	}
}

int main (void)
{
    std::vector<std::vector<double> > in;
    std::vector<std::vector<double> > R;
	std::vector<std::vector<double> > R_hated;
	std::vector<std::vector<int> > Ra;

    in = LoadInputMatrix();
    PrintMatrix(in, "Source matrix:\n");
    R = HemmingMatrix(in);
    PrintMatrix(R, "\nMatrix R:\n");

	R_hated = MinMaxTransitiveClosure(R);
	PrintMatrix(R_hated, "\nMatrix transitive closure of R:\n");

	RelationNegation(R_hated);

	Ra = AlphaLevel(R_hated, 1);



	getchar();
    return 0;
}
