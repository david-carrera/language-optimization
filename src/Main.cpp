#include <ParameterCollection.hpp>
#include <LexicalMatrixCSVGenerator.hpp>

int main(int argc, char *argv[]) {
    ParameterCollection params;
    params.parse(argc, argv);

    LexicalMatrixCSVGenerator csv(params);
    csv.generate();
    
    return EXIT_SUCCESS;
}
