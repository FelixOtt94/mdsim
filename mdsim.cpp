#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>
#include <errno.h>
#include <string.h>
#include <typeinfo>

// template class for parameters
class ParameterReader {

    private:
        std::map<std::string,std::string> data;

    public:

    bool readParameters(const std::string& filename){
        std::ifstream parameters (filename);
        std::string line;
        if (parameters.is_open()){
            while(!parameters.eof()){
                std::getline(parameters, line);
                if((line.length() <= 1)) continue;
                size_t pos = line.find_first_of(" ");
                size_t last = line.find_last_of(" ");
                data.insert(std::pair<std::string,std::string>(line.substr(0, pos), line.substr(last+1, line.length())));
            }
        }
        parameters.close();
        return true;
    }

    inline bool IsDefined(const std::string& key) const{
        int c = data.count(key);
        if(c == 0){
            return false;
        }
        return true;
    }

    template<typename Type>
    inline void GetParameter(const std::string& key, Type &value) const{
        std::map<std::string,std::string>::iterator it = data.find(key);
        if(typeid(value) == typeid(double)){
            value = std::stod(it->second);
        }
        else if(typeid(value) == typeid(int)){
            value = std::stoi(it->second);
        }
        else if(typeid(value) == typeid(std::string)){
            value = (it->second);
        }
    }

};

typedef struct {
   double m;
   double x0;
   double x1;
   double x2;
   double v0;
   double v1;
   double v2;
   double force0;
   double force1;
   double force2;
   double force0_old;
   double force1_old;
   double force2_old;
} particle;

void readInitialData(const char* filename, particle** particles, int* counter){
    double m, x0, x1, x2, v0, v1, v2;
    std::ifstream infile (filename);
    std::string line;
    if (infile.is_open())
    {
        while(!infile.eof()){
            std::getline(infile, line);
            if(!(line.length() <= 1))
                (*counter)++;
        }
    }
    infile.close();

    std::ifstream input (filename);

    (*particles) = (particle*)calloc(*counter, sizeof(particle));
    if(*particles == NULL){
        perror("malloc particles");
        exit(EXIT_FAILURE);
    }

    //memset((void*)particles, 0, *counter * sizeof(particle));

    if (input.is_open())
    {
        for(int i=0; i<*counter; ++i){
            input >> m >> x0 >> x1 >> x2 >> v0 >> v1 >> v2;

            (*particles)[i].m = m;
            (*particles)[i].x0 = x0;
            (*particles)[i].x1 = x1;
            (*particles)[i].x2 = x2;
            (*particles)[i].v0 = v0;
            (*particles)[i].v1 = v1;
            (*particles)[i].v2 = v2;
            //std::cout << m << " " << x0 << " " << x1 << " " << x2 << " " << v0 << " " << v1 << " " << v2 << std::endl;
        }
    }
    input.close();
}

void writeVTK(particle* particles, std::string &base, int time, int numParticles){
    std::string name(base);
    name += time;
    name += ".vtk";
    std::ofstream output (name);
    if(output.is_open())
    {
        output << "# vtk DataFile Version 3.0" << std::endl;
        output << "SiWiRVisFile" << std::endl;
        output << "ASCII" << std::endl;
        output << "DATASET UNSTRUCTURED_GRID" << std::endl;
        output << "POINTS " << numParticles << " DOUBLE" << std::endl;
        for(int i=0; i<numParticles; ++i){
            output << particles[i].x0 << " " << particles[i].x1 << " " << particles[i].x2 << std::endl;
        }
        output << "POINT_DATA " << numParticles << std::endl;
        output << "SCALARS mass double" << std::endl;
        output << "LOOKUP_TABLE default" << std::endl;
        for(int i=0; i<numParticles; ++i){
            output << particles[i].m << std::endl;
        }
        output << "VECTORS force double" << std::endl;
        for(int i=0; i<numParticles; ++i){
            output << particles[i].force0 << " " << particles[i].force1 << " " << particles[i].force2 << std::endl;
        }
        output << "VECTORS velocity double" << std::endl;
        for(int i=0; i<numParticles; ++i){
            output << particles[i].v0 << " " << particles[i].v1 << " " << particles[i].v2 << std::endl;
        }
    }
}

void simulation(particle* particles, int numParticles, ParameterReader &parameters){

    double t = 0, t_end = 0, delta_t = 0;
    int vis_space = 0;
    std::string name;
    parameters.GetParameter<double>(std::string("t_start"), t);
    parameters.GetParameter<double>(std::string("t_end"), t_end);
    parameters.GetParameter<double>(std::string("delta_t"), delta_t);
    parameters.GetParameter<int>(std::string("vis_space"), vis_space);
    parameters.GetParameter<std::string>(std::string("name"), name);
    int counter = 0;

    // compute forces Fi

    while(t < t_end){
        t += delta_t;
        for(int i=0; i<numParticles; ++i){
            particles[i].x0 += particles[i].v0 * delta_t + 0.5 * (particles[i].force0 / particles[i].m) * delta_t * delta_t;
            particles[i].x1 += particles[i].v1 * delta_t + 0.5 * (particles[i].force1 / particles[i].m) * delta_t * delta_t;
            particles[i].x2 += particles[i].v2 * delta_t + 0.5 * (particles[i].force2 / particles[i].m) * delta_t * delta_t;
            particles[i].force0_old = particles[i].force0;
            particles[i].force1_old = particles[i].force1;
            particles[i].force2_old = particles[i].force2;
        }

        // compute new forces Fi

        for(int i=0; i<numParticles; ++i){
            particles[i].v0 += particles[i].v0 + ((particles[i].force0_old + particles[i].force0)/(2*particles[i].m)) * delta_t;
            particles[i].v1 += particles[i].v1 + ((particles[i].force1_old + particles[i].force1)/(2*particles[i].m)) * delta_t;
            particles[i].v2 += particles[i].v2 + ((particles[i].force2_old + particles[i].force2)/(2*particles[i].m)) * delta_t;
        }
        counter++;
        if(counter%vis_space == 0){
            writeVTK(particles, name, counter, numParticles);
        }
    }

}

int main( int argc, char** argv ){
    if( argc != 3 ){
        std::cout << "Usage: ./mdsim [parameter file] [data file]" << std::endl;
        exit( EXIT_SUCCESS );
    }

    int numParticles = 0;
    particle* particles;

    readInitialData(argv[1], &particles, &numParticles);

    ParameterReader parameters;

    parameters.readParameters(std::string(argv[2]));

    simulation(particles, numParticles, parameters);

    free(particles);

    exit( EXIT_SUCCESS );
}
