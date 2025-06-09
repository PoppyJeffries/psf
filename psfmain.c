/***************************************************************************************
 * Filename: psfmain.c
 *
 * Description: Contains main entry-point into PSF simulator.
 ***************************************************************************************/

#include "psf.h"

/***************************************************************************************
 * Name: main()
 *
 * Purpose: Main function into PSF (Personality in Seabird Foraging).
 *
 * Parameters: (see usage statement below)
 *
 * Returns: zero
 ***************************************************************************************/
int main(int argc, char** argv)
{
	/* Local variables. */	
	int curr_arg;
	double total_time;
    unsigned short no_indivs;
    FILE *env_qual_file;
    FILE *patch_land_file;
    FILE *lambert_w_file;
    FILE *individual_output_file;
    FILE *patch_output_file;
    FILE *fs_energy_file;
    FILE *final_output_file;
    FILE *input_file;
    unsigned short no_lambert_w_inputs;
    double box_width;
    double box_height;
    double patch_app_rate;
    double mean_patch_longevity;
    double mean_patch_quality;
    double mean_step_length;
    double sensing_radius;
    double start_energy_level;
    double flying_energy_loss;
    double timestep;
    double cog_map_mult;
    double cog_map_decay;
    double cog_dir_mult;
    unsigned long timestamp;
    PSF_ROOT_DATA *root_data;
    
    /* If no arguments then print usage statement */
	if(argc == 1)
	{
		printf("Usage: psf.exe [-i <number-of-individuals> (default = 1)]\n");
		printf("               [-tt <total-time> (default = 100000)]\n");
		printf("               [-eq <environmental-quality-file>]\n");
		printf("               [-plf <patch-landing-prob-file>]\n");
		printf("               [-lwf <lambert-w-file>]\n");
		printf("               [-iof <individual_output_file>]\n");
    	printf("               [-pof <patch_output_file>]\n");
    	printf("               [-fsef <fs_energy_file>]\n");
    	printf("               [-fof <final_output_file>]\n");
    	printf("               [-if <input_file>]\n");
		printf("               [-nlw <no-lambert-w-inputs>]\n");
		printf("               [-bw <box-width> (default = 100)]\n");
		printf("               [-bh <box-height> (default = 100)]\n");
		printf("               [-par <patch-appearance-rate> (default = 0.1)]\n");
		printf("               [-mpl <mean-patch-longevity> (default = 10)]\n");
		printf("               [-mpq <mean-patch-quality> (default = 10)]\n");
		printf("               [-msl <mean-step-length> (default = 5)]\n");
		printf("               [-sr <sensing-radius> (default = 5)]\n");
		printf("               [-sel <start-energy-level> (default = 10)]\n");
		printf("               [-elf <energy-loss-flying> (default = 1)]\n");
		printf("               [-ts <timestep> (default = 0.1)]\n");
		printf("               [-cmm <cognitive-map-multiplier> (default = 0.1)]\n");
		printf("               [-cmd <cognitive-map-decay> (default = 0.1)]\n");
		printf("               [-cdm <cognitive-direction-multiplier> (default = 0.1)]\n");
		printf("               [-tst <timestamp> (default = unix time)]\n");
		goto EXIT_LABEL;
	}
	
	/* Default values */
    total_time = 100000;
    no_indivs = 1;
    box_width = 100;
    box_height = 100;
    no_lambert_w_inputs = 490;
    patch_app_rate = 0.1;
    mean_patch_longevity = 10;
    mean_patch_quality = 10;
    mean_step_length = 5;
    sensing_radius = 5;
    start_energy_level = 10;
    flying_energy_loss = 1;
    timestep = 0.1;
    cog_map_mult = 0.1;
    cog_map_decay = 0.1;
    cog_dir_mult = 0.1;
    
    timestamp = (unsigned long)time(NULL);

	
	/* Process arguments */
    for(curr_arg = 1; curr_arg < argc; curr_arg++)
    {
        if(!strcmp(argv[curr_arg], "-i"))
        {
            /* No of individuals */                      
         	no_indivs = atoi(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-tt"))
        {
            /* Total time */
        	total_time = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-bw"))
        {
            /* Box width */
        	box_width = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-bh"))
        {
            /* Box height */
        	box_height = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-par"))
        {
            /* Patch appearance rate */
        	patch_app_rate = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-mpl"))
        {
            /* Mean patch longevity */
        	mean_patch_longevity = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-mpq"))
        {
            /* Mean patch quality */
        	mean_patch_quality = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-msl"))
        {
            /* Mean step length */
        	mean_step_length = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-sr"))
        {
            /* Sensing radius */
        	sensing_radius = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-sel"))
        {
            /* Individual's starting energy level */
        	start_energy_level = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-elf"))
        {
            /* Flight energy loss */
        	flying_energy_loss = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-ts"))
        {
            /* Time between successive updates of the system */
        	timestep = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-tst"))
        {
            /* Timestamp */
        	timestamp = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-cmm"))
        {
            /* Cognitive map multiplier */
        	cog_map_mult = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-cmd"))
        {
            /* Cognitive map decay */
        	cog_map_decay = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-cdm"))
        {
            /* Cognitive direction multiplier */
        	cog_dir_mult = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-eq"))
        {
            /* Environment quality file */
        	env_qual_file = fopen(argv[++curr_arg],"r");
        }
        else if(!strcmp(argv[curr_arg], "-plf"))
        {
            /* File containing patch landing thresholds */
        	patch_land_file = fopen(argv[++curr_arg],"r");
        }
        else if(!strcmp(argv[curr_arg], "-lwf"))
        {
            /* File containing Lambert W function values */
        	lambert_w_file = fopen(argv[++curr_arg],"r");
        }
        else if(!strcmp(argv[curr_arg], "-nlw"))
        {
            /* Number of lambert w inputs */
        	no_lambert_w_inputs = atof(argv[++curr_arg]);
        }
        else if(!strcmp(argv[curr_arg], "-iof"))
        {
            /* File for individual output */
        	individual_output_file = fopen(argv[++curr_arg],"w");
        } 
		else if(!strcmp(argv[curr_arg], "-pof"))
        {
            /* File for patch output */
        	patch_output_file = fopen(argv[++curr_arg],"w");
        }
        else if(!strcmp(argv[curr_arg], "-fsef"))
        {
            /* File for patch output */
        	fs_energy_file = fopen(argv[++curr_arg],"w");
        }
		else if(!strcmp(argv[curr_arg], "-fof"))
        {
            /* File for patch output */
        	final_output_file = fopen(argv[++curr_arg],"a");
        }  
		else if(!strcmp(argv[curr_arg], "-if"))
        {
            /* File for patch output */
        	input_file = fopen(argv[++curr_arg],"a");
        }             
        else
        {
            /* Unrecognised parameter */
            fprintf(stderr,"Error: unrecognised parameter: %s\n", argv[curr_arg]);
            goto EXIT_LABEL;
        }
	}

	/* Initialise random seed. */
	srand(time(NULL)); 
	//srand(0); // Remove this after testing - keeps random numbers the same every time

	/* Initialisation: Enter values into CBs and allocate memory where necessary */
    root_data = psf_init(no_indivs,
                         total_time,
                         env_qual_file,
                         patch_land_file,
                         lambert_w_file,
                         individual_output_file,
                         patch_output_file,
                         fs_energy_file,
                         input_file,
                         final_output_file,
                         no_lambert_w_inputs,
		  				 box_height,
						 box_width,
						 patch_app_rate,
						 mean_patch_longevity,
						 mean_patch_quality,
						 mean_step_length,
						 sensing_radius,
						 start_energy_level,
						 flying_energy_loss,
						 cog_map_mult,
						 cog_map_decay,
						 cog_dir_mult,
						 timestep,
						 timestamp);
	if(root_data == NULL)
	{
		printf("Memory allocation failure\n");
		goto EXIT_LABEL;
	}

// START: Output for testing
    printf("Initialised\n");
// END: Output for testing
	
	/* Close the environmental quality file */
	if(env_qual_file != NULL)
	{
		fclose(env_qual_file);
	}

// START: Output for testing
    printf("Env file closed\n");
// END: Output for testing

	/* Close the patch land threshold file */
	if(patch_land_file != NULL)
	{
		fclose(patch_land_file);
	}

// START: Output for testing
    printf("Patch land file closed\n");
// END: Output for testing

	/* Close the patch land threshold file */
	if(lambert_w_file != NULL)
	{
		fclose(lambert_w_file);
	}

// START: Output for testing
    printf("Lambert W function file closed\n");
// END: Output for testing

// START: Output for data collection
    fprintf(root_data->input_file, "%u\t Input stats: \t -tt\t%f\t -bh\t%f\t -bw\t%f\t -par\t%f\t -mpl\t%f\t -mpq\t%f\t  -msl\t%f\t -sr\t%f\t -sel\t%f\t -elf\t%f\t -cmm\t%f\t -cmd\t%f\t -cmd\t%f\t -ts\t%f\t -i \t%i\n",
										root_data->timestamp,     
										root_data->total_time,
    									root_data->box_height,
										root_data->box_width,
										root_data->patch_app_rate,
										root_data->mean_patch_longevity,
										root_data->mean_patch_quality,
										root_data->mean_step_length,
										root_data->sensing_radius,
										root_data->start_energy_level,
										root_data->flying_energy_loss,
										root_data->cog_map_mult,
									    root_data->cog_map_decay,
									    root_data->cog_dir_mult,
										root_data->timestep,
										root_data->no_of_individuals
										);
// END: Output for data collection

	/* Run simulation */		
    psf_perform_simulation(root_data);
	
// START: Output for testing
    printf("Simulation done\n");
// END: Output for testing

// START: Output for testing
//    printf("Final individual data:\n");
// END: Output for testing

	/* Free memory allocated */
	psf_term(root_data);

// START: Output for testing
    printf("Termination complete\n");
// END: Output for testing

	/* Close the individual output file */
	if(individual_output_file != NULL)
	{
		fclose(individual_output_file);
	}
	
	/* Close the patch output file */
	if(patch_output_file != NULL)
	{
		fclose(patch_output_file);
	}

	/* Close the patch output file */
	if(fs_energy_file != NULL)
	{
		fclose(fs_energy_file);
	}

EXIT_LABEL:
	TRACE_END
	return(0);
}


