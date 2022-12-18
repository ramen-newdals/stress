double calcWSS_Parallel(int numPoints, int numThreads)
{
    double start_time, end_time;
    int i, j, k;
    double *wallShearRate = (double*)malloc(sizeof(double)*3);
    double *wallShearRateArray = (double*)malloc(sizeof(double)*(3*numPoints));
    double *shearVector = (double*)malloc(sizeof(double)*3);
    double *strainRateTensor = (double*)malloc(sizeof(double)*9);
    double *normal = (double*)malloc(sizeof(double)*3);
    double normalShear;
    //initialize Strain rate tensor
    for(i = 0; i<9; i++)
    {strainRateTensor[i] = i*i;}
    //Initialize norrmal vector
    for(i = 0; i<3; i++)
    {normal[i]=0.33;}
    omp_set_num_threads(numThreads);
    start_time = getElapsedTime();

}

double calcStrainRate_Parallel(int numPoints, int numThreads)
{
    double start_time, end_time;
    double *strainRateTensor = (double*)malloc(sizeof(double)*9);
    double *strainRateTensorArray = (double*)malloc(sizeof(double)*9*numPoints);
    double *velocityGradient = (double*)malloc(sizeof(double)*9);
    for(int i = 0; i<9; i++)
    {
        velocityGradient[i] = i*i;
    }
    start_time = getElapsedTime();
    omp_set_num_threads(numThreads);
    #pragma omp parallel
    for(int i=0; i<numPoints; i++)
    {
        for (int j=0; j<3; j++)
        {
            for (int k=0; k<3; k++)
            {
                strainRateTensor[3*j + k] = 0.5 * (velocityGradient[3*j + k] + velocityGradient[3*k + j]);
            }
        }
        for(int k = 0; k<9; k++)
        {
            strainRateTensorArray[(i*9)+k] = strainRateTensor[k];
        }
    }
    end_time = getElapsedTime();
    return (end_time - start_time);
}