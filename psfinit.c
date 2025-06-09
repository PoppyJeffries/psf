/***************************************************************************************
 * Filename: psfinit.c
 *
 * Description: Initialisation and termination procedures for the PSF simulator.
 ***************************************************************************************/

#include "psf.h"

/***************************************************************************************
 * Name: psf_init()
 *
 * Purpose: Allocate memory for CBs required for PSF simulator.  Set up initial CB
 *          constants.
 *
 * Parameters: IN     no_indivs - number of individuals in simulation.
 *             IN     total_time - total time for the simulation.
 *             IN     env_qual_file - file containing environmental quality.
 *             IN     patch_land_file - file containing details of patch landing thresholds.
 *             IN     lambert_w_file - file containing input and output values for the lambert w function.
 *             IN     individual_output_file - file to print individual output to
 *             IN     patch_output_file - file to print patch output to
 *             IN     fs_energy_file - file to print forage search energy outputs to
 *             IN     input_file - file to print input details to
 *             IN     final_output_file - file to print final output to
 *             IN     box_height - height of box where simulation is performed.
 *             IN     box_width - width of box where simulation is performed.
 *             IN     patch_app_rate - patch appearance rate.
 *             IN     mean_patch_longevity - mean patch longevity.
 *             IN     mean_patch_quality - mean patch quality.
 *             IN     mean_step_length - mean step length.
 *             IN     sensing_radius - sensing radius for patches.
 *             IN     start_energy_level - initial energy level of each bird.
 *             IN     flying_energy_loss - additional energy lost per unit time when flying.
 *             IN     cog_map_mult - cognitive map is updated by patch_quality*cog_map_multiplier.
 *             IN     cog_map_decay - cognitive map decay rate.
 *             IN     cog_dir_mult - cognitive direction multiplier in step selection weight function.
 *             IN     timestep - length of time step.
 *
 * Returns: root_data - pointer to a CB containing all the root data for PSF.  If
 *                      any memory allocation fail then return NULL.
 *
 * Operation: Allocate memory for root_data.  Input initial values.  For each of
 *            no_indivs individuals, call into the initialisation function for that
 *            individual.
 ***************************************************************************************/
PSF_ROOT_DATA *psf_init(unsigned short no_indivs,
                        double total_time,
                        FILE *env_qual_file,
                        FILE *patch_land_file,
                        FILE *lambert_w_file,
                        FILE *individual_output_file,
                        FILE *patch_output_file,
                        FILE *fs_energy_file,
                        FILE *input_file,
                        FILE *final_output_file,
                        unsigned short no_lambert_w_inputs,
		  				double box_height,
						double box_width,
						double patch_app_rate,
						double mean_patch_longevity,
						double mean_patch_quality,
						double mean_step_length,
						double sensing_radius,
						double start_energy_level,
						double flying_energy_loss,
						double cog_map_mult,
						double cog_map_decay,
						double cog_dir_mult,
						double timestep,
						unsigned long timestamp)	
{
    /* Local variables */
    unsigned short curr_individual;
    PSF_INDIVIDUAL *curr_indiv_cb;
    PSF_ROOT_DATA *root_data;
    unsigned long x_val;
    unsigned long y_val;
	
	/* Sanity checks. */
	assert(env_qual_file != NULL);
	
    /* Allocate memory for root data */
    root_data = (PSF_ROOT_DATA *)malloc(sizeof(PSF_ROOT_DATA));
	if(root_data == NULL)
	{
		fprintf(stderr,"MEMORY ALLOCATION FAILURE FOR ROOT_DATA\n");
		goto EXIT_LABEL;
	}

    /* Initial conditions of all the variables stored on root */
    root_data->total_time = total_time;
    root_data->box_height = box_height;
	root_data->box_width = box_width;
	root_data->patch_app_rate = patch_app_rate;
	root_data->mean_patch_longevity = mean_patch_longevity;
	root_data->mean_patch_quality = mean_patch_quality;
	root_data->mean_step_length = mean_step_length;
	root_data->sensing_radius = sensing_radius;
	root_data->start_energy_level = start_energy_level;
	root_data->flying_energy_loss = flying_energy_loss;	
	root_data->cog_map_mult = cog_map_mult;
    root_data->cog_map_decay = cog_map_decay;
    root_data->cog_dir_mult = cog_dir_mult;
    root_data->no_lambert_w_inputs = no_lambert_w_inputs;
	root_data->timestep = timestep;	
	root_data->no_of_individuals = no_indivs;
	root_data->patch_appearance_time = 0;
	root_data->current_time = 0;
	root_data->individual_output_file = individual_output_file;
	root_data->patch_output_file = patch_output_file;
	root_data->fs_energy_file = fs_energy_file;
	root_data->final_output_file = final_output_file;
	root_data->input_file = input_file;
	root_data->timestamp = timestamp;
// START: Output for testing
    printf("root_data->total_time:\t%f\n",root_data->total_time);
    printf("root_data->current_time:\t%f\n",root_data->current_time);
    printf("root_data->box_height:\t%f\n",root_data->box_height);
    printf("root_data->box_width:\t%f\n",root_data->box_width);
    printf("root_data->patch_app_rate:\t%f\n",root_data->patch_app_rate);
    printf("root_data->mean_patch_longevity:\t%f\n",root_data->mean_patch_longevity);
    printf("root_data->mean_patch_quality:\t%f\n",root_data->mean_patch_quality);
    printf("root_data->mean_step_length:\t%f\n",root_data->mean_step_length);
    printf("root_data->sensing_radius:\t%f\n",root_data->sensing_radius);
    printf("root_data->start_energy_level:\t%f\n",root_data->start_energy_level);
    printf("root_data->flying_energy_loss:\t%f\n",root_data->flying_energy_loss);
    printf("root_data->cog_map_mult:\t%f\n",root_data->cog_map_mult);
    printf("root_data->cog_map_decay:\t%f\n",root_data->cog_map_decay);
    printf("root_data->cog_dir_mult:\t%f\n", root_data->cog_dir_mult);
    printf("root_data->no_lambert_w_inputs:\t%f\n",root_data->no_lambert_w_inputs);
    printf("root_data->timestep:\t%f\n",root_data->timestep);
    printf("root_data->patch_appearance_time:\t%f\n",root_data->patch_appearance_time);
    printf("root_data->no_of_individuals:\t%i\n",root_data->no_of_individuals);
    printf("root_data->timestamp:\t%u\n",root_data->timestamp);
// END: Output for testing
    
    /* Get memory for patch landing thresholds */ 
	root_data->patch_land_thresh = (double *)calloc(no_indivs, sizeof(double));
	if(root_data->patch_land_thresh == NULL)
	{
		fprintf(stderr,"MEMORY ALLOCATION FAILURE FOR ROOT_DATA->patch_land_thresh\n");
		psf_term(root_data);
		root_data = NULL;
		goto EXIT_LABEL;
	}
	/* Put patch landing thresholds into array */
	if(patch_land_file != NULL)
	{
        psf_get_vector_from_file(root_data->patch_land_thresh, patch_land_file, no_indivs);
    }

    /* Initialise individuals */
    for(curr_individual = 0; curr_individual < no_indivs; curr_individual++)
    {
    	curr_indiv_cb = psf_indiv_init(curr_individual, root_data);
// START: Output for testing
        printf("Individual initialised at location %i\n", curr_indiv_cb);
// END: Output for testing
    	if(curr_indiv_cb == NULL)
    	{
            fprintf(stderr,"Failed to initialise individual %u\n", curr_individual);
            psf_term(root_data);
            root_data = NULL;
            goto EXIT_LABEL;
        }
    }

// START: Output for testing
    printf("Individuals initialised\n");
// END: Output for testing

    /* Get memory for environmental data */
	root_data->env_qual = (double *)calloc((unsigned long)(root_data->box_width*root_data->box_height), sizeof(double));
	if(root_data->env_qual == NULL)
	{
		fprintf(stderr,"MEMORY ALLOCATION FAILURE FOR ROOT_DATA->env_qual\n");
		psf_term(root_data);
		root_data = NULL;
		goto EXIT_LABEL;
	}
	/* Put environmental data into array */
	if(env_qual_file != NULL)
	{
        psf_get_array_from_file(root_data->env_qual, env_qual_file, root_data);
    }
    root_data->tot_env_qual = 0;
    for(y_val = 0; y_val < root_data->box_height; y_val++)
    {
        for(x_val = 0; x_val < root_data->box_width; x_val++)
        {
        	root_data->tot_env_qual += root_data->env_qual[x_val+y_val*(unsigned long)root_data->box_width];
// START: Output for testing
			//printf("Env quality at\t%i\t%i\tis\t%f\n", x_val, y_val, root_data->env_qual[x_val+y_val*(unsigned long)root_data->box_width]);
            //printf("Total environmental quality is\t%f\n", root_data->tot_env_qual);
// END: Output for testing
        }
    } 
    if(root_data->tot_env_qual == 0)
    {
		fprintf(stderr,"Cannot have an environment whose mean quality is zero\n");
		psf_term(root_data);
		root_data = NULL;
		goto EXIT_LABEL;
	}

// START: Output for testing
    printf("Environmental data set up\n");
    printf("Total environmental quality is\t%f\n", root_data->tot_env_qual);
// END: Output for testing
	    
	/* Get memory for Lambert W function */
	root_data->lambert_w = (double *)calloc((unsigned long)(root_data->no_lambert_w_inputs*2), sizeof(double));
	if(root_data->lambert_w == NULL)
	{
		fprintf(stderr,"MEMORY ALLOCATION FAILURE FOR ROOT_DATA->lambert_w\n");
		psf_term(root_data);
		root_data = NULL;
		goto EXIT_LABEL;
	}
	/* Put environmental data into array */
	if(lambert_w_file != NULL)
	{
        psf_get_vector_from_file(root_data->lambert_w, lambert_w_file, root_data->no_lambert_w_inputs*2);
    }

// START: Output for testing
    printf("Lambert W function array set up\n");
// END: Output for testing    
	    
EXIT_LABEL:
    /* Return pointer to the root data */
    return(root_data);
}

/***************************************************************************************
 * Name: psf_indiv_init()
 *
 * Purpose: Allocate memory for an individual.  Set up initial CB values.
 *
 * Parameters: IN     index - index for the individual
 *             IN     root_data - pointer to the per-PSF-simulator data.
 *
 * Returns: individual - pointer to a CB containing all the per-individual data.
 *
 * Operation: Allocate memory for individual.  Add the individual to the individual-list
 *            rooted in root_data.  Input initial values.
 ***************************************************************************************/
PSF_INDIVIDUAL *psf_indiv_init(unsigned short index, PSF_ROOT_DATA *root_data)
{
    /* Local variables */
    PSF_INDIVIDUAL *individual;
	unsigned long y_val;
	unsigned long x_val;

	/* Sanity checks.  The index must be between 0 and the number of individuals in the
	 * simulation.
	 */
	assert(root_data != NULL);
// START: Output for testing
    printf("index:\t%i\n",index);
// END: Output for testing
	assert(index < root_data->no_of_individuals);

	/* Allocate memory for the individual CB */
	individual = (PSF_INDIVIDUAL *)malloc(sizeof(PSF_INDIVIDUAL));
// START: Output for testing
    printf("Individual allocated at %i\n",individual);
// END: Output for testing
	if(individual == NULL)
	{
		fprintf(stderr, "MEMORY ALLOCATION FAILURE FOR INDIVIDUAL\n");
		goto EXIT_LABEL;
	}

	/* Add individual to list in root_data */
    LIST_V2_ADD_TO_START(&root_data->individual_list_root,
         	             &individual->list_elt,
    	                 individual);

// START: Output for testing
    printf("Individual added to list\n");
// END: Output for testing

	/* Input initial values */
    individual->current_x_pos = ((double)rand())*root_data->box_width/((double)RAND_MAX);
    individual->current_y_pos = ((double)rand())*root_data->box_height/((double)RAND_MAX);
    individual->initial_x = individual->current_x_pos;
    individual->initial_y = individual->current_y_pos;
    individual->index = index;
    individual->energy = root_data->start_energy_level;
    individual->patch_land_thresh = root_data->patch_land_thresh[index];
    individual->state = PSF_SEARCHING;
    individual->patch = NULL;       // Individual is not in a patch
    individual->ind_patch_count = 0;
    individual->energy_gain_count = 0;
    individual->energy_loss_count = 0;
    individual->search_start_time = 0;
    individual->avg_search_time = 0;

// START: Output for testing
    printf("Individual initial values assigned\n");
    printf("individual->current_x_pos:\t%f\n",individual->current_x_pos);
    printf("individual->current_y_pos:\t%f\n",individual->current_y_pos);
    printf("individual->index:\t%i\n",individual->index);
    printf("individual->energy:\t%f\n",individual->energy);
    printf("individual->patch_land_thresh:\t%f\n",individual->patch_land_thresh);
    printf("individual->state:\t%i\n",individual->state);
    printf("individual->patch:\t%i\n",individual->patch);
    printf("individual->ind_patch_count:\t%i\n",individual->ind_patch_count);
    printf("individual->energy_gain_count:\t%i\n",individual->energy_gain_count);
    printf("individual->energy_loss_count:\t%i\n",individual->energy_loss_count);
    printf("individual->search_start_time:\t%f\n",individual->search_start_time);
    printf("individual->avg_search_time:\t%f\n",individual->avg_search_time);
// END: Output for testing

    /* Set up cognitive map as array of zeros */	
    individual->cog_map = (double *)calloc((unsigned long)(root_data->box_width*root_data->box_height), sizeof(double));
// START: Output for testing
    printf("Cognitive map memory allocated at %i\n", individual->cog_map);
// END: Output for testing
	if(individual->cog_map == NULL)
	{
		fprintf(stderr,"MEMORY ALLOCATION FAILURE FOR ROOT_DATA->env_qual\n");
		psf_indiv_term(individual);
		individual = NULL;
		goto EXIT_LABEL;
	}
    for(y_val = 0; y_val < root_data->box_height; y_val++)
    {
        for(x_val = 0; x_val < root_data->box_width; x_val++)
        {
        	individual->cog_map[x_val+y_val*(unsigned long)(root_data->box_width)] = 0;
        }
    }     

EXIT_LABEL:
	/* Return pointer to the individual */
	return(individual);
}

/***************************************************************************************
 * Name: psf_patch_init()
 *
 * Purpose: Allocate memory for a patch.  Set up initial CB values.
 *
 * Parameters: IN     root_data - pointer to the per-PSF-simulator data.
 *
 * Returns: patch - pointer to a CB containing all the per-patch data.
 *
 * Operation: Allocate memory for patch.  Add the patch to the patch-list rooted in 
 *            root_data.  Input initial values.
 ***************************************************************************************/
PSF_PATCH *psf_patch_init(PSF_ROOT_DATA *root_data)
{
    /* Local variables */
    PSF_PATCH *patch;
    unsigned long x_val;
    unsigned long y_val;
    double random_no;
    double env_sum;
    unsigned short found_pos;

	/* Sanity checks. */
	assert(root_data != NULL);

	/* Allocate memory for the patch CB */
	patch = (PSF_PATCH *)malloc(sizeof(PSF_PATCH));
	if(patch == NULL)
	{
		fprintf(stderr, "MEMORY ALLOCATION FAILURE FOR PATCH\n");
		goto EXIT_LABEL;
	}

	/* Add individual to list in root_data */
    LIST_V2_ADD_TO_START(&root_data->patch_list_root,
         	             &patch->list_elt,
    	                 patch);

    /* Patch position is determined by underlying environment */
    random_no = ((double)rand())*root_data->tot_env_qual/((double)RAND_MAX + 1);
    env_sum = 0;
    found_pos = 0;
    for(y_val = 0; y_val < root_data->box_height; y_val++)
    {
        for(x_val = 0; x_val < root_data->box_width; x_val++)
        {
        	env_sum += root_data->env_qual[x_val+y_val*(unsigned long)(root_data->box_width)];
        	if(env_sum > random_no)
        	{
               patch->x_pos = x_val;
               patch->y_pos = y_val;
               found_pos = 1;
               break;
			}
        }
        if(found_pos == 1)
        {
        	break;
		}
    }

	/* Patch quality is drawn from an exponential distribuion */
	random_no = ((double)rand())/((double)RAND_MAX + 1);
    patch->quality = (0-log(1-random_no))*root_data->mean_patch_quality;
    
	/* Longevity is drawn from an exponential distribuion */
	random_no = ((double)rand())/((double)RAND_MAX + 1);
// START: Output for testing
//    printf("random_no: %f\n", random_no);
// END: Output for testing
    patch->longevity = (0-log(1-random_no))*root_data->mean_patch_longevity;
    
    /* Patch age is zero initially */
    patch->age = 0;
    patch->index = root_data->patch_count;

EXIT_LABEL:
	/* Return pointer to the individual */
	return(patch);
}
  
/***************************************************************************************
 * Name: psf_get_array_from_file()
 *
 * Purpose: Put values from a file into a box_width*box_height array.
 *
 * Parameters: OUT    array - array to copy values into
 *             IN     file  - file containing values
 *             IN     root_data - root data
 *
 * Returns: PSF_RC_OK/ERROR
 *
 * Operation: Read numbers from file and copy into array.  Numbers are stored as doubles.
 ***************************************************************************************/
void psf_get_array_from_file(double *array, FILE *file, PSF_ROOT_DATA *root_data)
{
	/* Local variables */
	char curr_char;
	char prev_char;
	double curr_val = 0;
	unsigned long y_val = 0;
	unsigned long x_val = 0;
	unsigned short got_to_minus = PSF_NO;
	double digit_after_dec = 0;

	/* Sanity checks */
	assert(array != NULL);
    assert(file != NULL);

    rewind(file);
	curr_char = getc(file);
	while(curr_char != EOF)
	{
		if((curr_char >= '0') && (curr_char <= '9') && (digit_after_dec == 0) && (got_to_minus == PSF_NO))
        { 
            /* This is a digit of a number */
            curr_val = curr_val*10+(double)(curr_char-'0');
        }
  		else if((curr_char >= '0') && (curr_char <= '9') && (got_to_minus == PSF_NO))
        { 
            /* This is a digit of a number after the decimal point */
            curr_val += (double)(curr_char-'0')/pow(10,digit_after_dec);
            digit_after_dec++;
        }
        else if((curr_char == '.') && (got_to_minus == PSF_NO))
        {
            /* Decimal point. */
            digit_after_dec++;
        }
        else if(curr_char == '-')
        {
             /* Minus sign.  Number is 0 */
             curr_val = 0;
             got_to_minus = PSF_YES;
        }
		else if((curr_char == '\n') && (prev_char >= '0') && (prev_char <= '9'))
		{
            /* End of row of numbers */
            assert(x_val <= root_data->box_width);
            assert(y_val <= root_data->box_height);
            array[x_val+y_val*(unsigned long)(root_data->box_width)] = curr_val;
            curr_val = 0;
            y_val++;
            x_val = 0;
            digit_after_dec = 0;
            got_to_minus = PSF_NO;
        } 
		else if((curr_char == '\n') && (prev_char == '\t'))
		{
            /* End of row of numbers after a tab.  No need to store values in array.  Just
             * change x_val and y_val to reflect the fact that we are about to start a new
             * row.
             */
            assert(x_val+y_val*root_data->box_width <= root_data->box_height*root_data->box_width);
            curr_val = 0;
            y_val++;
            x_val = 0;
            digit_after_dec = 0;
        } 
    	else if((curr_char == '\t') && (prev_char >= '0') && (prev_char <= '9'))
		{
            /* End of a number */
            assert(x_val <= root_data->box_width);
            assert(y_val <= root_data->box_height);
            array[x_val+y_val*(unsigned long)(root_data->box_width)] = curr_val;
            curr_val = 0;
            x_val++;
            got_to_minus = PSF_NO;
            digit_after_dec = 0;
        }                            
        else
        {
            /* Unrecognised character.  Maybe benign.  No-op */
        }             
        prev_char = curr_char;              
    	curr_char = getc(file);
	}
	
	/* Return */
	return;
}

/***************************************************************************************
 * Name: psf_get_vector_from_file()
 *
 * Purpose: Put values from a file into a 1D array of given length.
 *
 * Parameters: OUT    array - array to copy values into
 *             IN     file  - file containing values
 *             IN     length - check if all values in array are less than this.  
 *                                nb. could be PSF_INFINITY
 *
 * Returns: Nothing
 *
 * Operation: Read numbers from file and copy into array.  Numbers are stored as doubles.
 ***************************************************************************************/
void psf_get_vector_from_file(double *array, FILE *file, unsigned long length)
{
	/* Local variables */
	char curr_char;
	char prev_char;
	double curr_val = 0;
	unsigned long x_val = 0;
	double digit_after_dec = 0;
	signed char plus_minus = 1;

	/* Sanity checks */
	assert(array != NULL);
    assert(file != NULL);

    rewind(file);
	curr_char = getc(file);
	while(curr_char != EOF)
	{
		if((curr_char >= '0') && (curr_char <= '9') && (digit_after_dec == 0))
        { 
            /* This is a digit of a number */
            curr_val = curr_val*10+(double)(curr_char-'0');
        }
  		else if((curr_char >= '0') && (curr_char <= '9'))
        { 
            /* This is a digit of a number after the decimal point */
            curr_val += (double)(curr_char-'0')/pow(10,digit_after_dec);
            digit_after_dec++;
        }
        else if(curr_char == '.')
        {
            /* Decimal point. */
            digit_after_dec++;
        }
        else if(curr_char == '-')
        {
             /* Minus sign. */
             plus_minus = -1;
        }
		else if((curr_char == '\n') && (prev_char >= '0') && (prev_char <= '9'))
		{
            /* End of row of numbers */
            curr_val *= plus_minus;
            assert(x_val <= length);
            array[x_val] = curr_val;
// START: Output for testing
            //printf("array value: %f\n", curr_val);
// END: Output for testing
            curr_val = 0;
            x_val++;
            digit_after_dec = 0;
            plus_minus = 1;
        } 
    	else if((curr_char == '\t') && (prev_char >= '0') && (prev_char <= '9'))
		{
            /* End of a number */
            curr_val *= plus_minus;
            assert(x_val <= length);
            array[x_val] = curr_val;
// START: Output for testing
            //printf("array value: %f\n", curr_val);
// END: Output for testing
            curr_val = 0;
            x_val++;
            digit_after_dec = 0;
            plus_minus = 1;
        }                            
        else
        {
            /* Unrecognised character.  Maybe benign.  No-op */
        }             
        prev_char = curr_char;              
    	curr_char = getc(file);
    	
	}
	
	/* Return */
	return;
}

/***************************************************************************************
 * Name: psf_term()
 *
 * Purpose: Free all memory allocated in psf_init.
 *
 * Parameters: IN     root_data - pointer to a CB containing all the root data for PSF.
 *
 * Returns: Nothing.
 ***************************************************************************************/
void psf_term(PSF_ROOT_DATA *root_data)
{
    /* Local variables */
    PSF_INDIVIDUAL *curr_indiv_cb;
    PSF_PATCH *curr_patch_cb;

	/* Sanity checks. */
	assert(root_data != NULL);

    /* Free arrays */
    if(root_data->env_qual != NULL)
    {
        free(root_data->env_qual);
	}
    if(root_data->patch_land_thresh != NULL)
    {
        free(root_data->patch_land_thresh);
	}
    
    /* Free up patch CBs */
    if(!LIST_V2_EMPTY(root_data->patch_list_root))
    {
        curr_patch_cb = (PSF_PATCH *)LIST_V2_GET_FIRST(root_data->patch_list_root);
        while(curr_patch_cb != NULL)
        {
            /* Delete the first element on the list.  Then free the associated control block.
             * Then get the new first element. */
            LIST_V2_DELETE_FIRST(root_data->patch_list_root);
    	    psf_patch_term(curr_patch_cb);
    	    curr_patch_cb = (PSF_PATCH *)LIST_V2_GET_FIRST(root_data->patch_list_root);
        }
	}

    /* Free up individual CBs */
    if(!LIST_V2_EMPTY(root_data->individual_list_root))
    {
        curr_indiv_cb = (PSF_INDIVIDUAL *)LIST_V2_GET_FIRST(root_data->individual_list_root);
        while(curr_indiv_cb != NULL)
        {
            /* Delete the first element on the list.  Then free the associated control block.
             * Then get the new first element. */
            LIST_V2_DELETE_FIRST(root_data->individual_list_root);
    	    psf_indiv_term(curr_indiv_cb);
    	    curr_indiv_cb = (PSF_INDIVIDUAL *)LIST_V2_GET_FIRST(root_data->individual_list_root);
        }
    }
    
    /* Free root */
	free(root_data);
        
EXIT_LABEL:
    /* Return */
    return;
}

/***************************************************************************************
 * Name: psf_indiv_term()
 *
 * Purpose: Free all memory allocated in psf_indiv_init().
 *
 * Parameters: IN     individual - pointer to per-individual data.
 *
 * Returns: Nothing.
 *
 * Operation: Free per-individual memory.
 ***************************************************************************************/
void psf_indiv_term(PSF_INDIVIDUAL *individual)
{
	/* Local variables */

	/* Sanity checks */
	assert(individual != NULL);

	/* Free up the memory allocated for various arrays */
	free(individual->cog_map);

	/* Free up the memory allocated for the individual control block */
	free(individual);

	/* Return */
	return;
}

/***************************************************************************************
 * Name: psf_patch_term()
 *
 * Purpose: Free all memory allocated in psf_patch_init().
 *
 * Parameters: IN     patch - pointer to per-patch data.
 *
 * Returns: Nothing.
 *
 * Operation: Free per-individual memory.
 ***************************************************************************************/
void psf_patch_term(PSF_PATCH *patch)
{
	/* Local variables */

	/* Sanity checks */
	assert(patch != NULL);

	/* Free up the memory allocated for the patch control block */
	free(patch);

	/* Return */
	return;
}
