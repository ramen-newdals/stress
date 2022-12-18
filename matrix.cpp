void constructLSMatricies(int matrixSize)
{
    Eigen::MatrixXd A(matrixSize, matrixSize);
    Eigen::VectorXd b(matrixSize);
    Eigen::VectorXd x(matrixSize);
    double avg_length=0;
    int i, j;
    for(i = 0; i<matrixSize; i++)
    {
        b(i) = i;
        for(j = 0; j<matrixSize; j++)
        {
            A(i, j) = i*j;
        }
    }
    x = A.colPivHouseholderQr().solve(b);
    //std::cout << "Hello From: " << omp_get_thread_num() << std::endl;
}


double constructMatriciesParallel(int matrixSize)
{
    double start_time, end_time, run_time;
    double *A = (double*)malloc(sizeof(double)*(matrixSize*matrixSize));
    start_time = getElapsedTime();
    #pragma omp parallel
    for(int i = 0; i<(matrixSize*matrixSize); i++)
    {
        A[i] = i*i;
    }
    end_time = getElapsedTime();
    run_time = (end_time - start_time);
    return run_time;
}

double constructMatriciesSerial(int matrixSize)
{
    double start_time, end_time, run_time;
    Eigen::MatrixXd A(matrixSize, matrixSize);
    start_time = getElapsedTime();
    for(int i = 0; i<matrixSize; i++)
    {
        for(int j = 0; j<matrixSize; j++)
        {
            A(i, j) = i*j;
        }
    }
    end_time = getElapsedTime();
    run_time = (end_time - start_time);
    return run_time;
}