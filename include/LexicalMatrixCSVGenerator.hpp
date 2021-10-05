#ifndef LEXICAL_MATRIX_CSV_GENERATOR_H
#define LEXICAL_MATRIX_CSV_GENERATOR_H

#include <Util.hpp>
#include <LexicalMatrixParameters.hpp>
#include <LexicalMatrixStatistics.hpp>
#include <BaseLexicalMatrix.hpp>
#include <LexicalMatrixOptimizer.hpp>
#include <ParameterCollection.hpp>
#include <DataCollector.hpp>

/**
 * Generate CSVs of a LexicalMatrix, calling the optimizer, etc.
 */
class LexicalMatrixCSVGenerator {
private:
    /// Header used for the information theoretical CSV file
    const std::string info_header = "lambda\tH(S)\tH(R)\tH(S,R)\tI(S,R)\tH(S|R)\tH(R|S)\tRefUseless\tConnectComp\tConnectCompWords\tConnectCompMeanings\n";
    /// Header used for the lambda specific statistical information CSV file
    const std::string lambda_header = "p_i\tmu_i\ta_i\tmu\n";

    /// For updating prcess status, last time the status was updated as to avoid flooding the screen too fast.
    std::chrono::time_point<std::chrono::steady_clock> lastStatusUpdate;
    /// Time when the process started
    std::chrono::time_point<std::chrono::steady_clock> startTime;
    /// Value of the lambda being currently processed
    double currentLambda;
    /// Number of the realization being currently processed
    unsigned currentSample;

    /// Parameters to use for the process
    ParameterCollection &params;

    /// Data collectors used to generate CSVs.
    std::vector<BaseDataCollector*> data_collectors;

    /// Show current status of the process to the user
    void update_status();

    /// Set data collectors to 0
    void reset_collectors();

    void collect_collectors(const BaseLexicalMatrix &mat,
                            const LexicalMatrixOptimizer &opti);
    void record_collectors(double lambda);
public:
    LexicalMatrixCSVGenerator(ParameterCollection &params);
    ~LexicalMatrixCSVGenerator();
    void generate();
};

#endif /* LEXICAL_MATRIX_CSV_GENERATOR_H */
