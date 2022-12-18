double readMesh(int numTrials)
{
    double start_time, end_time, avg_time=0;
    for(int i = 0; i<numTrials; i++)
    {
        meshReader mesh1;
        start_time = getElapsedTime();
        mesh1.readMesh();
        end_time = getElapsedTime();
        avg_time += (end_time - start_time);
    }
    return (avg_time/numTrials);
}