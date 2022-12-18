std::vector<float> compareParallelSerial()
{
    double start_time, end_time;
    int i;
    std::vector<float> solutionTimes; 
    start_time = getElapsedTime();
    omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions
    #pragma omp parallel
    for(i = 0; i<8; i++)
    {constructLSMatricies(10);}
    end_time = getElapsedTime();
    std::cout << "Parallel Sollution Took: " << end_time - start_time << std::endl;
    solutionTimes.push_back((end_time-start_time));
    start_time = getElapsedTime();
    for(int i=0; i<8; i++)
    {constructLSMatricies(10);}
    end_time = getElapsedTime();
    std::cout << "Serial Sollution Took: " << end_time - start_time << std::endl;
    solutionTimes.push_back((end_time-start_time));
    return solutionTimes;
}
