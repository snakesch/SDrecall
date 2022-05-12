#include <iostream>
#include <cstdlib>
#include <thread>
#include <cstring>

using namespace std;

int set(void)
{
    /* Optimal number of threads for computation */
    unsigned int maxThreads = thread::hardware_concurrency();
    char* reqThreads = to_string(maxThreads / 4);
    setenv("_THREADS", reqThreads, true);
}

char* int2charptr(int i)
{
    char *ret = NULL;
    while (i > 0)
    {
        char cur = (char) (i % 10 + '0');
        i /= 10;
    }
}

int unset(void)
{
    unsetenv("_THREADS");
}

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cerr << "Error: No instructions given!\n";
        return -1;
    }
    if ( strcmp(*(argv + 1), "set") == 0 )
        set();
    else if (strcmp(*(argv + 1), "unset") == 0)
        unset();
    else
    {
        cerr << "Error: Unknown instruction " << *(argv + 1) << endl;
        return -1;
    }

    return 0;
}