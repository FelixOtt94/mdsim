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
#include <list>
#include <vector>

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

	  
	  
     inline void GetParameter(const std::string& key, std::string &value) const{
            std::map<std::string,std::string>::const_iterator it = data.find(key);
	      value = (it->second);
     }
     
         template<typename Type>
    inline void GetParameter(const std::string& key, Type &value) const{
        std::map<std::string,std::string>::const_iterator it = data.find(key);
        if(typeid(value) == typeid(double)){
            value = std::stod(it->second);
	    return;
        }
        else if(typeid(value) == typeid(int)){
            value = std::stoi(it->second);
	    return;
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

typedef struct{
    double size_cell_x;
    double size_cell_y;
    double size_cell_z;
    int numbers_cell_x;
    int numbers_cell_y;
    int numbers_cell_z;
}cell;

void readInitialData(const char* filename, particle** particles, int* counter, std::vector<std::list<int>> & cell_array, ParameterReader& parameters, cell* cell_parameter){
    double m, x0, x1, x2, v0, v1, v2;
    double size_x, size_y, size_z;
    double x_min, y_min, z_min;
    double r_cut;
    int x_index, y_index, z_index;

    parameters.GetParameter<double>(std::string("x_max"), size_x);
    parameters.GetParameter<double>(std::string("x_min"), x_min);
    size_x -= x_min;
    parameters.GetParameter<double>(std::string("y_max"), size_y);
    parameters.GetParameter<double>(std::string("y_min"), y_min);
    size_y -= y_min;
    parameters.GetParameter<double>(std::string("z_max"), size_z);
    parameters.GetParameter<double>(std::string("z_min"), z_min);
    size_z -= z_min;

    parameters.GetParameter<double>(std::string("r_cut"), r_cut);

    cell_parameter-> numbers_cell_x = ((int) (size_x / r_cut));
    cell_parameter-> size_cell_x = size_x/((double) cell_parameter-> numbers_cell_x);
    cell_parameter-> numbers_cell_y = ((int) (size_y / r_cut));
    cell_parameter-> size_cell_y = size_y/((double) cell_parameter-> numbers_cell_y);
    cell_parameter-> numbers_cell_z = ((int) (size_z / r_cut));
    cell_parameter-> size_cell_z = size_z/((double) cell_parameter-> numbers_cell_z);

    cell_array.resize((cell_parameter->numbers_cell_x) * (cell_parameter->numbers_cell_y) * (cell_parameter->numbers_cell_z));


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

            x_index = (int)((x0-x_min)/cell_parameter -> size_cell_x);
            y_index = (int)((x1-y_min)/cell_parameter -> size_cell_y);
            z_index = (int)((x2-z_min)/cell_parameter -> size_cell_z);

            cell_array[z_index*(cell_parameter->numbers_cell_z) + y_index*(cell_parameter->numbers_cell_y) + x_index].push_front(i);
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
    parameters.GetParameter<double>( std::string("t_start"), t);
    parameters.GetParameter<double>( std::string("t_end"), t_end);
    parameters.GetParameter<double>( std::string("delta_t"), delta_t);
    parameters.GetParameter<int>( std::string("vis_space"), vis_space);
    parameters.GetParameter( std::string("name"), name);
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

void updateForce(ParameterReader &parameters, particle *particles, cell *cell_parameter, std::vector<std::list<int>> &cell_array){

    int neighbors [27];
    double r_ij;
    double r [3];
    double r_cut, epsilon, sigma;
    parameters.GetParameter<double>(std::string("r_cut"), r_cut);
    parameters.GetParameter<double>(std::string("epsilon"), epsilon);
    parameters.GetParameter<double>(std::string("sigma"), sgima);

    const double epsilon_24 = 24.0*epsilon;
    const double sigma_6 = sigma*sigma*sigma*sigma*sigma*sigma;

    for(std::vector<std::list>::iterator it_ic = cell_array.begin(); it_ic != cell_array.end(); ++it_ic){
        for(std::list::iterator it_i = cell_array[*it_ic].begin(); it_i != cell_array[*it_ic].end(); ++it_i){
            particles[*it_i].force0 = 0.0;
            particles[*it_i].force1 = 0.0;
            particles[*it_i].force2 = 0.0;

            //neighbors()
            for(int i = 0; i<27; ++i){
                for(std::list::iterator it_j = cell_array[neighbors[i]].begin(); it_j != cell_array[neighbors[i]].end(); ++it_j){
                        if(*it_i != *it_j){
                            r[0]=particles[*it_j].x0 - particles[*it_i].x0;
                            r[1]=particles[*it_j].x1 - particles[*it_i].x1;
                            r[2]=particles[*it_j].x2 - particles[*it_i].x2;
                            r_ij = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
                            if(r_ij <= r_cut){
                                particles[*it_i].force0 += epsilon_24 * (1.0/(r_ij*r_ij))XXXXXXXXXXXXXX;
                                particles[*it_i].force1 += ;
                                particles[*it_i].force2 += ;
                            }
                        }

                }
            }

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
    cell cell_parameter;
    std::vector<std::list<int>> cell_array;


    ParameterReader parameters;

    parameters.readParameters(std::string(argv[2]));

    readInitialData(argv[1], &particles, &numParticles, cell_array, parameters, &cell_parameter);



    simulation(particles, numParticles, parameters);

    free(particles);

    exit( EXIT_SUCCESS );
}
