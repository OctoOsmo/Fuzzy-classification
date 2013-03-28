#include <stdio.h>
#include <vector>
#include <cmath>

std::vector<std::vector<double> > LoadInputMatrix(/*std::vector<std::vector<double> > &in*/)
{
    std::vector<std::vector<double> > in;
    //Чтение из файла исходной информации
    int size = 0;
    FILE *f = fopen("in", "r+");
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
    return l;
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

int main (void)
{
    std::vector<std::vector<double> > in;
    std::vector<std::vector<double> > R;
    in = LoadInputMatrix();
    PrintMatrix(in, "Исходная матрица:\n");
    R = HemmingMatrix(in);
    PrintMatrix(R, "Матрица R:\n");
    return 0;
}
