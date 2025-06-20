/***************************************************************************************
 * Filename: psf.h
 *
 * Description: Header file containing all PSF CBs, defines, includes, macros and
 *              function declarations that are required in PSF.
 ***************************************************************************************/

#ifndef __PSF_H_
#define __PSF_H_

/***************************************************************************************
 * Includes
 ***************************************************************************************/
#include <assert.h>
#include <stdio.h>
#include <time.h> 
#include <stdlib.h>
//#include <listjrp.h>
#include "listv2.h" /* flag this up */
#include <math.h>
#include <stdbool.h>

/***************************************************************************************
 * Pre-processor definitions
 ***************************************************************************************/

/***************************************************************************************
 * Yes and no.
 ***************************************************************************************/
#define PSF_YES 1
#define PSF_NO  2

/***************************************************************************************
 * Behavioural state: searching or foraging
 ***************************************************************************************/
#define PSF_SEARCHING     1
#define PSF_FORAGING      2

/***************************************************************************************
 * Number of samples for determining next move
 ***************************************************************************************/
#define PSF_NO_SAMPLES     100

/***************************************************************************************
 * Approx value of pi
 ***************************************************************************************/
#define PSF_PI (3.14159265358979)

/***************************************************************************************
 * Macros
 ***************************************************************************************/

/***************************************************************************************
 * Calculate distance from a line to a point
 *
 * Parameters: X0 - x-value of the point
 *             Y0 - y-value of the point
 *             X1 - x-value of point one on the line
 *             Y1 - y-value of point one on the line
 *             X2 - x-value of point two on the line
 *             Y2 - y-value of point two on the line
 ***************************************************************************************/
#define PSF_DIST_LINE_TO_POINT(X0,Y0,X1,Y1,X2,Y2) \
    fabs(((Y2)-(Y1))*(X0)-((X2)-(X1))*(Y0)+(X2)*(Y1)-(X1)*(Y2))/sqrt(pow((Y2)-(Y1),2)+pow((X2)-(X1),2)) 
//#define PSF_DIST_LINE_TO_POINT(X0,Y0,X1,Y1,X2,Y2) \
    ((Y2)-(Y1))*(X0)-((X2)-(X1))*(Y0)+(X2)*(Y1)-(X1)*(Y2)/(sqrt(pow((Y2)-(Y1),2)+pow((X2)-(X1),2)))

/***************************************************************************************
 * Control blocks
 ***************************************************************************************/

/***************************************************************************************
 * Name: PSF_PATCH
 *
 * Purpose: Stores all data relating to resources.
 ***************************************************************************************/
typedef struct psf_patch
{
	/***********************************************************************************
	 * x-position
	 ***********************************************************************************/
    double x_pos;

	/***********************************************************************************
	 * y-position.
	 ***********************************************************************************/
    double y_pos;

	/***********************************************************************************
	 * List element of this patch.
	 ***********************************************************************************/
    LIST_V2_ELT list_elt;

	/***********************************************************************************
	 * Patch quality
	 ***********************************************************************************/
    double quality;
    
	/***********************************************************************************
	 * Patch longevity
	 ***********************************************************************************/
    double longevity;

	/***********************************************************************************
	 * Patch age
	 ***********************************************************************************/
    double age;
    
    /***********************************************************************************
	 * Patch index
	 ***********************************************************************************/
    unsigned long index;

} PSF_PATCH;

/***************************************************************************************
 * Name: PSF_INDIVIDUAL
 *
 * Purpose: Stores all individual-level data.
 ***************************************************************************************/
typedef struct psf_individual
{
	/***********************************************************************************
	 * Current x-position
	 ***********************************************************************************/
    double current_x_pos;

	/***********************************************************************************
	 * Current y-position.
	 ***********************************************************************************/
    double current_y_pos;

	/***********************************************************************************
	 * List element of this individual.
	 ***********************************************************************************/
    LIST_V2_ELT list_elt;

	/***********************************************************************************
	 * Index of this individual.
	 ***********************************************************************************/
    unsigned short index;

	/***********************************************************************************
	 * Patch that the individual is in (NULL if searching).
	 ***********************************************************************************/
    PSF_PATCH *patch;
    
    /***********************************************************************************
	 * Number of patches visited by this individual.
	 ***********************************************************************************/
    unsigned long ind_patch_count;
    
    /***********************************************************************************
	 * Number of times energy is gained in a forage-search event.
	 ***********************************************************************************/
    unsigned long energy_gain_count;
    
    /***********************************************************************************
	 * Number of times energy is lost in a forage-search event.
	 ***********************************************************************************/
    unsigned long energy_loss_count;
    
    /***********************************************************************************
	 * Searching start time, used to calculate avg_search_time.
	 ***********************************************************************************/
    double search_start_time;
    
    /***********************************************************************************
	 * Average search time of individual.
	 * Equivalent to 1/Lambda from Jeffries et al. model.
	 ***********************************************************************************/
    double avg_search_time;  

    /***********************************************************************************
	 * Foraging start time, used to calculate time to leave patch.
	 ***********************************************************************************/
    double foraging_start_time;
    
    /***********************************************************************************
	 * Time to leave patch.
	 ***********************************************************************************/
    double opt_foraging_time;
    
	/***********************************************************************************
	 * Average patch quality for individual.
	 ***********************************************************************************/
    double avg_patch_quality;
    
    /***********************************************************************************
	 * Initial quality of the patch the individual is currently foraging in.
	 ***********************************************************************************/
    double curr_patch_initial_quality;
    
	/***********************************************************************************
	 * Initial quality of the patch the individual is currently foraging in.
	 ***********************************************************************************/
    double curr_patch_leaving_thres;
    
    /***********************************************************************************
	 * Initial energy of individual when joining a patch
	 ***********************************************************************************/
    double patch_start_energy;
    
	/***********************************************************************************
	 * Details of previous patch for updating cognitive map after a new patch has been found
	 ***********************************************************************************/
    double prev_patch_x;
	double prev_patch_y;
	double prev_patch_quality;

	/***********************************************************************************
	 * Average straight line distance between patches
	 ***********************************************************************************/

    double avg_distance_between_patches;
    
    /***********************************************************************************
	 * Details for distance travelled
	 ***********************************************************************************/
    double initial_x;
    double initial_y;
    double distance_travelled_since_patch;
    double avg_distance_travelled_between_patches;

	/***********************************************************************************
	 * Energy level
	 ***********************************************************************************/
    double energy;

	/***********************************************************************************
	 * Patch landing threshold
	 ***********************************************************************************/
    double patch_land_thresh;

	/***********************************************************************************
	 * State: searching or foraging
	 ***********************************************************************************/
    unsigned short state;

	/***********************************************************************************
	 * Cognitive map
	 ***********************************************************************************/
    double *cog_map;
    
	/***********************************************************************************
	 * Target search location
	 ***********************************************************************************/
	 double target_x;
	 double target_y;
    
} PSF_INDIVIDUAL;

/***************************************************************************************
 * Name: PSF_ROOT_DATA
 *
 * Purpose: The place in which all PSF data is rooted.
 ***************************************************************************************/
typedef struct psf_root_data
{
	/***********************************************************************************
	 * Time for simulation to run.
	 ***********************************************************************************/
    double total_time;

	/***********************************************************************************
	 * Time between successive updates of the simulation.
	 ***********************************************************************************/
    double timestep;

	/***********************************************************************************
	 * Multiplier to apply to patch quality when updating the cognitive map (i.e. cog map
	 * is updated by cog_map_mult*patch->quality).
	 ***********************************************************************************/
    double cog_map_mult;
    
    /***********************************************************************************
	 * Decay rate of cognitive map
	 ***********************************************************************************/
    double cog_map_decay;
 
 	/***********************************************************************************
	 * Multiplier in step selection weight function for direction to target location
	 ***********************************************************************************/   
    double cog_dir_mult;

	/***********************************************************************************
	 * Number of individuals in simulation.
	 ***********************************************************************************/
    unsigned short no_of_individuals;

	/***********************************************************************************
	 * Root of the list of psf_individual CBs.
	 ***********************************************************************************/
    LIST_V2_ROOT individual_list_root;

	/***********************************************************************************
	 * Root of the list of psf_patch CBs.
	 ***********************************************************************************/
    LIST_V2_ROOT patch_list_root;

	/***********************************************************************************
	 * Current time.
	 ***********************************************************************************/
    double current_time;
    
	/***********************************************************************************
	 * Time last patch appeared.
	 ***********************************************************************************/
    double patch_appearance_time;
    
	/***********************************************************************************
	 * Patch appearance rate.
	 ***********************************************************************************/
    double patch_app_rate;
    
	/***********************************************************************************
	 * Mean patch longevity.
	 ***********************************************************************************/
    double mean_patch_longevity;
    
	/***********************************************************************************
	 * Mean patch quality.
	 ***********************************************************************************/
    double mean_patch_quality;

	/***********************************************************************************
	 * Time since last patch appearance.
	 ***********************************************************************************/
    double last_patch_time;
    
	/***********************************************************************************
	 * Mean step length.
	 ***********************************************************************************/
    double mean_step_length;

	/***********************************************************************************
	 * Sensing radius - for sensing the presence of a patch.
	 ***********************************************************************************/
    double sensing_radius;

	/***********************************************************************************
	 * Width of box.
	 ***********************************************************************************/
    double box_width;

	/***********************************************************************************
	 * Height of box.
	 ***********************************************************************************/
    double box_height;
    
	/***********************************************************************************
	 * Starting energy level for bird.
	 ***********************************************************************************/
    double start_energy_level;
    
	/***********************************************************************************
	 * Additional energy lost from flying per unit time.
	 ***********************************************************************************/
    double flying_energy_loss;
    
	/***********************************************************************************
	 * Array of patch landing thresholds
	 ***********************************************************************************/
    double *patch_land_thresh;
    
	/***********************************************************************************
	 * Array giving environmental quality
	 ***********************************************************************************/
    double *env_qual;
    
	/***********************************************************************************
	 * Array of Lambert W function values
	 ***********************************************************************************/
    double *lambert_w;
    
    /***********************************************************************************
	 * Number of Lambert W entries
	 ***********************************************************************************/
    unsigned short no_lambert_w_inputs;
    
	/***********************************************************************************
	 * Sum of all the env_qual values
	 ***********************************************************************************/
    double tot_env_qual;
    
	/***********************************************************************************
	 * File to write individual data to
	 ***********************************************************************************/    
    FILE *individual_output_file;

	/***********************************************************************************
	 * File to write patch data to
	 ***********************************************************************************/    
    FILE *patch_output_file;
	
	/***********************************************************************************
	 * File to write individual data to
	 ***********************************************************************************/    
    FILE *fs_energy_file;       
    
    /***********************************************************************************
	 * File to write input data to
	 ***********************************************************************************/    
    FILE *input_file;  
    
    /***********************************************************************************
	 * File to write final output data to
	 ***********************************************************************************/    
    FILE *final_output_file;  
    
    /***********************************************************************************
	 * Sum of all the env_qual values
	 ***********************************************************************************/
    unsigned long timestamp;
    
    /***********************************************************************************
	 * Patch counter
	 ***********************************************************************************/
    unsigned long patch_count;
    

} PSF_ROOT_DATA;

/***************************************************************************************
 * Function declarations.
 ***************************************************************************************/

/***************************************************************************************
 * psfmain.c
 ***************************************************************************************/
int main(int, char**);

/***************************************************************************************
 * psfinit.c
 ***************************************************************************************/
PSF_ROOT_DATA *psf_init(unsigned short,
                        double,
                        FILE *,
                        FILE *,
                        FILE *,
                        FILE *,
                        FILE *,
                        FILE *,
                        FILE *,
                        FILE *,                        
                        unsigned short,
						double,
						double,
						double,
						double,
						double,
						double,
						double,
						double,
						double,
						double,
						double,
						double,
						double,
						unsigned long);
PSF_INDIVIDUAL *psf_indiv_init(unsigned short, PSF_ROOT_DATA *);
PSF_PATCH *psf_patch_init(PSF_ROOT_DATA *);
void psf_get_array_from_file(double *, FILE *, PSF_ROOT_DATA *);
void psf_get_vector_from_file(double *, FILE *, unsigned long);
void psf_term(PSF_ROOT_DATA *);
void psf_indiv_term(PSF_INDIVIDUAL *);
void psf_patch_term(PSF_PATCH *);

/***************************************************************************************
 * psfproc.c
 ***************************************************************************************/
void psf_perform_simulation(PSF_ROOT_DATA *);
void psf_move_individual(PSF_INDIVIDUAL *, PSF_ROOT_DATA *);
double psf_dist_seg_to_point(double, 
                             double, 
							 double, 
							 double,
							 double, 
							 double);
double psf_lambert_w(double, PSF_ROOT_DATA *);
void psf_new_target_loc(PSF_INDIVIDUAL *, PSF_ROOT_DATA *);

#endif /* __PSF_H_ */
