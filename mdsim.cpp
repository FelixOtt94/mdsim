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
    double size_x;
    double size_y;
    double size_z;
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
    cell_parameter->size_x = size_x;
    cell_parameter->size_y = size_y;
    cell_parameter->size_z = size_z;

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
            //input >> m >> x0 >> x1 >> x2 >> v0 >> v1 >> v2;
             std::getline(input, line);
             int first = line.find_first_not_of(std::string(" "), 0);
             int seconde = line.find(" ", first );
             m = stod( line.substr(first, seconde) );
             first = line.find_first_not_of(std::string(" "), seconde);
             seconde = line.find(" ", first );
             x0 = stod( line.substr(first, seconde) );
             first = line.find_first_not_of(std::string(" "), seconde);
             seconde = line.find(" ", first );
             x1 = stod( line.substr(first, seconde) );
             first = line.find_first_not_of(std::string(" "), seconde);
             seconde = line.find(" ", first );
             x2 = stod( line.substr(first, seconde) );
             first = line.find_first_not_of(std::string(" "), seconde);
             seconde = line.find(" ", first );
             v0 = stod( line.substr(first, seconde) );
             first = line.find_first_not_of(std::string(" "), seconde);
             seconde = line.find(" ", first );
             v1 = stod( line.substr(first, seconde) );
             first = line.find_first_not_of(std::string(" "), seconde);
             seconde = line.find(" ", first );
             v2 = stod( line.substr(first, seconde) );

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

            cell_array[z_index*(cell_parameter->numbers_cell_x*cell_parameter->numbers_cell_y) + y_index*(cell_parameter->numbers_cell_x) + x_index].push_front(i);
        }
    }
    input.close();



}

void writeVTK(particle* particles, std::string &base, int time, int numParticles){
    std::string name(base);
    name += std::to_string(time);
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

void getNeighbors(int i, int* neighbors,  cell* cell_parameter){
    int counter =0;
    int zkod = i / (cell_parameter->numbers_cell_x * cell_parameter->numbers_cell_y);
    int ykod = (i - zkod*(cell_parameter->numbers_cell_x * cell_parameter->numbers_cell_y)) / cell_parameter->numbers_cell_x;
    int xkod = i - zkod*(cell_parameter->numbers_cell_x * cell_parameter->numbers_cell_y) - ykod*cell_parameter->numbers_cell_x;
    int xkod_up, ykod_up, zkod_up;
    //std::cout << "x: " << xkod << " y: " << ykod << " z: " << zkod << std::endl;
    for(int z=-1; z<2; z++){
        zkod_up = zkod+z;
        if(zkod+z < 0){
            zkod_up = cell_parameter->numbers_cell_z-1;
        }else if(zkod+z >= cell_parameter->numbers_cell_z ){
            zkod_up = 0;
        }
        for(int y=-1; y<2; y++){
            ykod_up = ykod+y;
            if(ykod+y < 0){
                ykod_up = cell_parameter->numbers_cell_y-1;
            }else if(ykod+y >= cell_parameter->numbers_cell_y ){
                ykod_up = 0;
            }
            for(int x=-1; x<2; x++){
                xkod_up = xkod+x;
                if(xkod+x < 0){
                    xkod_up = cell_parameter->numbers_cell_x-1;
                }else if(xkod+x >= cell_parameter->numbers_cell_x ){
                    xkod_up = 0;
                }
                //std::cout << "x: " << xkod_up << " y: " << ykod_up << " z: " << zkod_up << std::endl;
                //neighbors[counter++] = i + cell_parameter->numbers_cell_x *y + (cell_parameter->numbers_cell_x * cell_parameter->numbers_cell_y)*z+ x;
                neighbors[counter++] = xkod_up + ykod_up*cell_parameter->numbers_cell_x+ zkod_up*(cell_parameter->numbers_cell_x * cell_parameter->numbers_cell_y);
            }
        }
    }
}

void updateForce(ParameterReader &parameters, particle *particles, std::vector<std::list<int>> &cell_array, cell* cell_parameter){

    int neighbors [27];
    int c=-1;
    double r_ij;
    double r [3];
    double r_cut, epsilon, sigma;
    parameters.GetParameter<double>(std::string("r_cut"), r_cut);
    parameters.GetParameter<double>(std::string("epsilon"), epsilon);
    parameters.GetParameter<double>(std::string("sigma"), sigma);

    const double epsilon_24 = 24.0*epsilon;
    const double sigma_6 = sigma*sigma*sigma*sigma*sigma*sigma;

    for(std::vector<std::list<int>>::iterator it_ic = cell_array.begin(); it_ic != cell_array.end(); ++it_ic){
        c++;
        for(std::list<int>::iterator it_i = (*it_ic).begin(); it_i != (*it_ic).end(); ++it_i){
            particles[*it_i].force0 = 0.0;
            particles[*it_i].force1 = 0.0;
            particles[*it_i].force2 = 0.0;

            getNeighbors( c, neighbors , cell_parameter);
            for(int i = 0; i<27; ++i){
                for(std::list<int>::iterator it_j = cell_array[neighbors[i]].begin(); it_j != cell_array[neighbors[i]].end(); ++it_j){

                        if(*it_i != *it_j){

                            r[0]=particles[*it_j].x0 - particles[*it_i].x0;
                            r[1]=particles[*it_j].x1 - particles[*it_i].x1;
                            r[2]=particles[*it_j].x2 - particles[*it_i].x2;
                            if(r[0] > cell_parameter->size_cell_x){
                                r[0] =  -cell_parameter->size_x + fabs(particles[*it_j].x0 - particles[*it_i].x0);
                            }else if( r[0] < (-cell_parameter->size_cell_x) ){
                                r[0] =  cell_parameter->size_x - fabs(particles[*it_j].x0 - particles[*it_i].x0);
                            }
                            if(r[1] > cell_parameter->size_cell_y){
                                r[1] =  -cell_parameter->size_y + fabs(particles[*it_j].x1 - particles[*it_i].x1);
                            }else if( r[1] < (-cell_parameter->size_cell_y) ){
                                r[1] =  cell_parameter->size_y - fabs(particles[*it_j].x1 - particles[*it_i].x1);
                            }
                            if(r[2] > cell_parameter->size_cell_z){
                                r[2] =  -cell_parameter->size_z + fabs(particles[*it_j].x2 - particles[*it_i].x2);
                            }else if( r[2] < (-cell_parameter->size_cell_z) ){
                                r[2] =  cell_parameter->size_z - fabs(particles[*it_j].x2 - particles[*it_i].x2);
                            }
                            r_ij = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
                            if(r_ij <= r_cut){

                                double tmp = (epsilon_24 * (1.0/(r_ij*r_ij)) * (sigma_6/pow(r_ij,6)) * (1.0-2.0*(sigma_6/pow(r_ij,6))));
                                particles[*it_i].force0 += tmp * r[0];
                                particles[*it_i].force1 += tmp * r[1];
                                particles[*it_i].force2 += tmp * r[2];
                                //std::cout << particles[*it_i].force0 << " " << particles[*it_i].force1 << " " << std::endl;
                                if(*it_i == 8)
                                    std::cout << "rij " << r_ij << " r0 " << r[0] << " r1 " << r[1] << " tmp " << tmp << " f0 "<< particles[8].force0 <<  " f1 "<< particles[8].force1 << std::endl;
                            }
                        }

                }
            }

        }


    }

}

void simulation(particle* particles, int numParticles, ParameterReader &parameters, std::vector<std::list<int>> &cell_array, cell* cell_parameter){

    double t = 0, t_end = 0, delta_t = 0;
    int x_index, y_index, z_index;
    double x_min, y_min, z_min;
    int index_old, index_new;
    int vis_space = 0;
    std::string name;
    parameters.GetParameter<double>( std::string("t_start"), t);
    parameters.GetParameter<double>( std::string("t_end"), t_end);
    parameters.GetParameter<double>( std::string("delta_t"), delta_t);
    parameters.GetParameter<int>( std::string("vis_space"), vis_space);
    parameters.GetParameter( std::string("name"), name);
    parameters.GetParameter<double>(std::string("x_min"), x_min);
    parameters.GetParameter<double>(std::string("y_min"), y_min);
    parameters.GetParameter<double>(std::string("z_min"), z_min);
    int counter = 0;
std::cout << "PARTICLE 8: x: " << particles[8].x0 << " y: " << particles[8].x1 << " v0 " << particles[8].v0 << " v1 " << particles[8].v1 << " f0 " << particles[8].force0 << " f1 " << particles[8].force1 << std::endl;
    updateForce(parameters, particles, cell_array, cell_parameter);
std::cout << "PARTICLE 8: x: " << particles[8].x0 << " y: " << particles[8].x1 << " v0 " << particles[8].v0 << " v1 " << particles[8].v1 << " f0 " << particles[8].force0 << " f1 " << particles[8].force1 << std::endl;
    while(t < t_end){

        t += delta_t;
        for(int i=0; i<numParticles; ++i){
            x_index = (int)((particles[i].x0-x_min)/cell_parameter -> size_cell_x);
            y_index = (int)((particles[i].x1-y_min)/cell_parameter -> size_cell_y);
            z_index = (int)((particles[i].x2-z_min)/cell_parameter -> size_cell_z);
//std::cout << "koords old x: " << particles[i].x0 << " y: " << particles[i].x1 << " z: " << particles[i].x2 << std::endl;
//std::cout << "index old x: " << x_index << " y: " << y_index << " z: " << z_index << std::endl;
            index_old = z_index*(cell_parameter->numbers_cell_y*cell_parameter->numbers_cell_x) + y_index*(cell_parameter->numbers_cell_x) + x_index;
            if(i==8)
std::cout << "PARTICLE 8: x: " << particles[i].x0 << " y: " << particles[i].x1 << " v0 " << particles[i].v0 << " v1 " << particles[i].v1 << " f0 " << particles[i].force0 << " f1 " << particles[i].force1 << std::endl;
            particles[i].x0 += particles[i].v0 * delta_t + 0.5 * (particles[i].force0 / particles[i].m) * delta_t * delta_t;
            particles[i].x1 += particles[i].v1 * delta_t + 0.5 * (particles[i].force1 / particles[i].m) * delta_t * delta_t;
            particles[i].x2 += particles[i].v2 * delta_t + 0.5 * (particles[i].force2 / particles[i].m) * delta_t * delta_t;
            particles[i].force0_old = particles[i].force0;
            particles[i].force1_old = particles[i].force1;
            particles[i].force2_old = particles[i].force2;
            x_index = (int)((particles[i].x0-x_min)/cell_parameter -> size_cell_x);
            if( x_index >= cell_parameter->numbers_cell_x ){
                x_index = 0;
            }else if( x_index < 0 ){
                x_index = cell_parameter->numbers_cell_x-1;
            }
            y_index = (int)((particles[i].x1-y_min)/cell_parameter -> size_cell_y);
            if( y_index >= cell_parameter->numbers_cell_y ){
                y_index = 0;
            }else if( y_index < 0 ){
                y_index = cell_parameter->numbers_cell_y-1;
            }
            z_index = (int)((particles[i].x2-z_min)/cell_parameter -> size_cell_z);
            if( z_index >= cell_parameter->numbers_cell_z ){
                 z_index = 0;
            }else if( z_index < 0 ){
                z_index = cell_parameter->numbers_cell_z-1;
            }
//std::cout << "x: " << x_index << " y: " << y_index << " z: " << z_index << std::endl;
            index_new = z_index*(cell_parameter->numbers_cell_y*cell_parameter->numbers_cell_x) + y_index*(cell_parameter->numbers_cell_x) + x_index;
            if( index_old != index_new){
//std::cout << "index_old: " << index_old << " index_new " << index_new << std::endl;
                cell_array[index_old].remove(i);
                cell_array[index_new].push_front(i);
                if( particles[i].x0 >= (cell_parameter->size_x + x_min) ){
                    particles[i].x0 -= cell_parameter->size_x;
                }else if( particles[i].x0 <= x_min ){
                    particles[i].x0 += cell_parameter->size_x;
//if(particles[i].x0 < 0)
//std::cout << i << " "<< particles[i].v0 << std::endl;
                }
                if( particles[i].x1 >= (cell_parameter->size_y + y_min) ){
                    particles[i].x1 -= cell_parameter->size_y;
                }else if( particles[i].x1 <= y_min ){
                    particles[i].x1 += cell_parameter->size_y;
                }
                if( particles[i].x2 >= (cell_parameter->size_z + z_min) ){
                    particles[i].x2 -= cell_parameter->size_z;
                }else if( particles[i].x2 <= z_min ){
                    particles[i].x2 += cell_parameter->size_z;
                }
            }
        }
//std::cout << "upadte force" << std::endl;
        updateForce(parameters, particles, cell_array, cell_parameter);

        for(int i=0; i<numParticles; ++i){
            particles[i].v0 +=  ((particles[i].force0_old + particles[i].force0)/(2*particles[i].m)) * delta_t;
            particles[i].v1 +=  ((particles[i].force1_old + particles[i].force1)/(2*particles[i].m)) * delta_t;
            particles[i].v2 +=  ((particles[i].force2_old + particles[i].force2)/(2*particles[i].m)) * delta_t;
        }
        counter++;
        if(counter%vis_space == 0){
//std::cout << "vtk" << std::endl;
            writeVTK(particles, name, counter, numParticles);
        }
    }

}



int main( int argc, char** argv ){
    int n[27];
    cell p;
    p.numbers_cell_x = 3;
    p.numbers_cell_y = 3;
    p.numbers_cell_z = 3;
    getNeighbors(0, n,  &p);
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                std::cout << n[i*9+j*3+k] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

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


    simulation(particles, numParticles, parameters, cell_array, &cell_parameter);

    free(particles);

    exit( EXIT_SUCCESS );
}
