#include <LexicalMatrixCSVGenerator.hpp>
#include <LexicalMatrixStatistics.hpp>
#include <BaseLexicalMatrix.hpp>
#include <fstream>

void LexicalMatrixCSVGenerator::update_status() {
    if (!this->params.show_status()) {
        return;
    }
    
    auto currentTime = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedTime = currentTime - this->lastStatusUpdate;
    if (elapsedTime.count() < this->params.status_update_period()) {
        return;
    }
    this->lastStatusUpdate = currentTime;

    elapsedTime = currentTime - this->startTime;
    unsigned elapsedRuns = this->params.elapsed_runs(this->currentLambda) + this->currentSample;
    unsigned totalRuns = this->params.total_runs();
    unsigned remainingRuns = totalRuns - elapsedRuns;

    double timePerRun;
    if (elapsedRuns != 0) {
        timePerRun = elapsedTime.count() / elapsedRuns;
    } else {
        timePerRun = std::numeric_limits<double>::infinity();
    }
    
    double remainingTime = remainingRuns * timePerRun;
    SeparatedTime time(remainingTime);

    unsigned percentageCompleted;
    if (totalRuns != 0) {
        percentageCompleted = 100 * elapsedRuns / totalRuns;
    } else {
        percentageCompleted = 0;
    }
    
    std::cerr << " |l:" << this->params.elapsed_lambdas(this->currentLambda) << "/" << this->params.total_lambdas()-1;
    std::cerr << "|s:" << this->currentSample+1 << "/" << this->params.realizations();
    std::cerr << "| ";

    std::cerr << percentageCompleted << "% ";
    std::cerr << "ETA: " << time;

    std::cerr << "          ";
    
    std::cerr << std::flush << "\r";
}

void LexicalMatrixCSVGenerator::reset_collectors() {
    for (BaseDataCollector *collector : this->data_collectors) {
        collector->reset();
    }
}

/**
 * Extract data needed for data collectors from the matrix and the optimizer object.
 * @param[in]  mat   Lexical Matrix
 * @param[in]  opti  Optimizer
 */
void LexicalMatrixCSVGenerator::collect_collectors(const BaseLexicalMatrix &mat, const LexicalMatrixOptimizer &opti) {
    for (BaseDataCollector *collector : this->data_collectors) {
        collector->collect(mat, opti);
    }
}

/**
 * Save data from data collectors to disk.
 * @param[in]  lambda  Lambda used to generate the collected data
 */
void LexicalMatrixCSVGenerator::record_collectors(double lambda) {
    for (BaseDataCollector *collector : this->data_collectors) {
        collector->record(lambda);
    }
}

/**
 * Construct object from parameters.
 */
LexicalMatrixCSVGenerator::LexicalMatrixCSVGenerator(ParameterCollection &params) : params(params) {
    unsigned n = this->params.n();
    unsigned m = this->params.m();
    unsigned samples = this->params.realizations();
    
    if (this->params.generate_info()) {
        BaseDataCollector *collector = new InfoDataCollector();
        collector->setParams([&params](double lambda){
            (void)lambda;
            return params.info_csv();
        }, n, m, samples);
        this->data_collectors.push_back(collector);
    }

    if (this->params.generate_lambdas()) {
        BaseDataCollector *collector = new LambdaDataCollector();
        collector->setParams([&params](double lambda){
            return params.lambda_csv(lambda);
        }, n, m, samples);
        this->data_collectors.push_back(collector);
    }

    if (this->params.generate_graph()) {
        BaseDataCollector *collector = new GraphDataCollector();
        collector->setParams([&params](double lambda){
            return params.graphs_dir(lambda);
        }, n, m, samples);
        this->data_collectors.push_back(collector);
    }

    #ifdef TRACE_MATRIX_OPTIMIZATION
    {
        BaseDataCollector *collector = new OptimizationTraceDataCollector();
        collector->setParams([&params](double lambda){
            (void)lambda;
            return params.base_dir();
        }, n, m, samples);
        this->data_collectors.push_back(collector);
    }
    #endif
}

/**
 * Destroy object, removing all data collectors.
 */
LexicalMatrixCSVGenerator::~LexicalMatrixCSVGenerator() {
    for (BaseDataCollector *collector : this->data_collectors) {
        delete collector;
    }
}

/**
 * Start the process of generating all the data.
 */
void LexicalMatrixCSVGenerator::generate() {
    MultipleOstreams outInfoTxt;
    
    std::ofstream info_txt;
    info_txt.open(this->params.info_txt());
    if (!info_txt.fail()) {
        outInfoTxt.add(&info_txt);
    }
    
    if (this->params.show_status()) {
        outInfoTxt.add(&std::cerr);
    }

    outInfoTxt << "Parameters:\n";
    outInfoTxt << this->params << "\n";
    outInfoTxt.flushall();
    
    if (this->params.generate_lambdas()) {
        std::error_code ec;
        filesystem::remove_all(this->params.lambda_dir(), ec);
        filesystem::create_directories(this->params.lambda_dir());
    }

    if (this->params.generate_info()) {
        std::ofstream info_file;
        info_file.open(this->params.info_csv());
        if (info_file.fail()) {
            std::cerr << "Failed to open file: " << this->params.info_csv() << std::endl;
            return;
        }
        info_file << this->info_header;
        info_file.close();
    }

    LexicalMatrixOptimizer opti = this->params.make_optimizer();
    BaseLexicalMatrix &mat = this->params.make_matrix();

    this->startTime = std::chrono::steady_clock::now();
    
    auto startDate = std::chrono::system_clock::now();
    std::time_t startDateTime = std::chrono::system_clock::to_time_t(startDate);
    outInfoTxt << "Started on: " << std::ctime(&startDateTime);

    std::vector<double> lambdas;
    this->params.lambdas(lambdas);
    for (double lambda : lambdas) {
        this->currentLambda = lambda;

        if (this->params.generate_lambdas()) {
            std::ofstream lambdas_file;
            auto filename = this->params.lambda_csv(lambda);
            lambdas_file.open(filename);
            if (lambdas_file.fail()) {
                std::cerr << "Failed to open file: " << filename << std::endl;
                return;
            }
            lambdas_file << this->lambda_header;
            lambdas_file.close();
        }

        this->reset_collectors();
        
        for (unsigned j=0; j<this->params.realizations(); j++) {
            opti.optimize(mat, lambda);
            this->collect_collectors(mat, opti);
            this->currentSample = j;
            this->update_status();
        }

        if (this->params.generate_graph()) {
            auto graphs_dir = this->params.graphs_dir(lambda);
            std::error_code ec;
            filesystem::remove_all(graphs_dir, ec);
            filesystem::create_directories(graphs_dir);
        }
        
        this->record_collectors(lambda);
    }

    auto finishDate = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = finishDate - startDate;
    std::time_t finishDateTime = std::chrono::system_clock::to_time_t(finishDate);

    outInfoTxt << "Finished on: " << std::ctime(&finishDateTime);
    outInfoTxt << "Time elapsed: " << SeparatedTime(elapsed_seconds.count()) << "\n";
    outInfoTxt.flushall();
    
    info_txt.close();
}
