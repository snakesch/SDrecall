#include <zlib.h>
#include <iostream>
#include <cstring>
#include <cstdlib>

#define MAX_HEADER_LEN 400
#define MAX_READ_LEN 5000

using std::cout;
using std::cerr;
using std::endl;

void printError(const char* line1, const char* line2, const char* line3, const char* line4)
{
    cerr << "A snippet of the read is as follows: " << endl;
    cerr << line1
         << line2
         << line3
         << line4
         << endl;
    return;
}

bool isValidFQ(const char * path, unsigned long long expectedLn = 80000000)
{
    unsigned long long rd = 0;

    /* Open gz file */
    gzFile infile = gzopen(path, "rb");

    while (!gzeof(infile))
    {
        char line1[MAX_HEADER_LEN + 1] = {};
        char line2[MAX_READ_LEN + 1] = {};
        char line3[5] = {};
        char line4[MAX_READ_LEN + 1] = {};

        bool extraLine = gzeof(infile);
        gzgets(infile, line1, MAX_HEADER_LEN);
        if (gzeof(infile) && !extraLine)
            return 0;
        gzgets(infile, line2, MAX_READ_LEN);
        gzgets(infile, line3, 5);
        gzgets(infile, line4, MAX_READ_LEN);

        if (*line3 != '+')
            return 0;
        if (line1 == NULL || line2 == NULL || line3 == NULL || line4 == NULL)
            return 0;

        rd++;
    }

    /* Close gz file */
    gzclose(infile);
    
    return 1;
}

static void printUsage(char* argv0)
{
    cout << "Usage: " << argv0 << endl
                      << "\t-i inputPath -t fileType [-o outputPath] [-e expectedLineNo]" << endl;
    return;
}

int main(int argc, char** argv)
{
    int cur_pos = 1;
    char* inputFN = NULL;
    char* outputFN = NULL;
    char* fileType = NULL;
    unsigned long long expectedLn = 0;
    if (argc < 2)
    {
        printUsage(argv[0]);
        return -1;
    }
    while (cur_pos < argc)
    {
        if(strcmp(*(argv + cur_pos), "-h") == 0)
        {
            printUsage(argv[0]);
            return 0;
        }
        else if (strcmp(*(argv + cur_pos), "-i") == 0)
        {
            cur_pos++;
            inputFN = *(argv + cur_pos);
        }
        else if (strcmp(*(argv + cur_pos), "-o") == 0)
        {
            cur_pos++;
            outputFN = *(argv + cur_pos);
        }
        else if (strcmp(*(argv + cur_pos), "-e") == 0)
        {
            cur_pos++;
            expectedLn = strtoull(*(argv + cur_pos), NULL, 0);
        }
        else if (strcmp(*(argv + cur_pos), "-t") == 0)
        {
            cur_pos++;
            fileType = *(argv + cur_pos);
        }
        else
        {
            printUsage(argv[0]);
            return -1;
        }

        cur_pos++;
    }

    if (strcmp(fileType, "fastq") == 0 || strcmp(fileType, "fq") == 0)
    {
        /* Check file validity - emptiness / Ln = multiples of 4 */
        bool code = (expectedLn == 0) ? isValidFQ(inputFN) : isValidFQ(inputFN, expectedLn);
        cout << !code;
    }

    return 0;
}
