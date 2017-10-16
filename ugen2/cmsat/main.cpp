
#include <ctime>
#include <cstring>
#include <errno.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <map>
#include <set>
#include "cmsat/constants.h"
#include <signal.h>
#include <time.h>
#include <list>
#include <algorithm>
#include <signal.h>
#include <functional>

#include "cmsat/time_mem.h"
#include "cmsat/constants.h"
#include "cmsat/DimacsParser.h"

#if defined(__linux__)
#include <fpu_control.h>
#endif

#include "cmsat/Main.h"
#include "cmsat/main.h"

using namespace CMSat;

/**
@brief For correctly and gracefully exiting

It can happen that the user requests a dump of the learnt clauses. In this case,
the program must wait until it gets to a state where the learnt clauses are in
a correct state, then dump these and quit normally. This interrupt hander
is used to achieve this
 */

int main(int argc, char** argv) {
    Main main(argc, argv);
    main.parseCommandLine();
    signal(SIGINT, SIGINT_handler);
    //signal(SIGALRM, SIGALARM_handler);
    try{
        return main.singleThreadSolve();

    }

    catch(std::bad_alloc) {
        std::cerr << "Memory manager cannot handle the load. Sorry. Exiting." << std::endl;
        exit(-1);
    }

    catch(std::out_of_range oor) {
        std::cerr << oor.what() << std::endl;
        exit(-1);
    }

    catch(CMSat::DimacsParseError dpe) {
        std::cerr << "PARSE ERROR!" << dpe.what() << std::endl;
        exit(3);
    }
    return 0;
}
