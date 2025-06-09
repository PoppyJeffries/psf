/***************************************************************************************
 * Filename: psfproc.c
 *
 * Description: General procedures for the PSF simulator.
 ***************************************************************************************/

#include "psf.h"

/***************************************************************************************
 * Name: psf_perform_simulation()
 *
 * Purpose: Perform the simulation.
 *
 * Parameters: IN    root_data - pointer to a CB containing all the root data for PSF.
 *
 * Returns: Nothing.
 *
 * Operation: Go through individuals, moving each and updating their energy levels.  Add
 *            patches as a Poisson process.  Terminate patches if either longevity is 
 *            reached or energy is below 0.
 ***************************************************************************************/
void psf_perform_simulation(PSF_ROOT_DATA *root_data)
{
	/* Local variables */
	PSF_INDIVIDUAL *curr_indiv_cb;
	PSF_INDIVIDUAL *prev_indiv_cb;
	PSF_PATCH *curr_patch_cb;
	PSF_PATCH *prev_patch_cb;
	double random_no;
	double curr_cum_dist;
	unsigned long no_patches_added;
	unsigned long k_factorial;

    /* Sanity checks. Current time should be 0. */
    assert(root_data != NULL);
    assert(root_data->current_time == 0); 
	assert(root_data->patch_count ==0); 
	
	while(root_data->current_time <= root_data->total_time)
	{
        /* Increment time */
	    root_data->current_time += root_data->timestep;
// START: Output for testing
        //printf("Time:\t%f\n", root_data->current_time); 
// END: Output for testing

        /* Get first individual in list */
        curr_indiv_cb = (PSF_INDIVIDUAL *)LIST_V2_GET_FIRST(root_data->individual_list_root);
        /* If the list is empty then end the simulation */
        if(curr_indiv_cb == NULL)
        {
        	printf("All the birds have died! :(\n");
        	break;
		}
        /* Loop through individuals, moving each */
        while(curr_indiv_cb != NULL)
        {
        	/* Move this individual, updating its energy */
            psf_move_individual(curr_indiv_cb, root_data);
// START: Output for data collection
            fprintf(root_data->individual_output_file, "%u\t Individual %i:\t%f\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\n", 
			                                           root_data->timestamp,
													   curr_indiv_cb->index,
            										   curr_indiv_cb->patch_land_thresh,
													   root_data->current_time, 
			                                           curr_indiv_cb->current_x_pos, 
												       curr_indiv_cb->current_y_pos,
												       curr_indiv_cb->energy,
												       curr_indiv_cb->state,
													   curr_indiv_cb->ind_patch_count,
													   curr_indiv_cb->energy_gain_count,
													   curr_indiv_cb->energy_loss_count,
													   curr_indiv_cb->avg_patch_quality,
													   curr_indiv_cb->avg_search_time,
													   curr_indiv_cb->avg_distance_travelled_between_patches,
													   curr_indiv_cb->avg_distance_between_patches
													   );
// END: Output for data collection
			if(root_data->current_time > root_data->total_time)
			{
// START: Output for data collection
        		printf("Individual %i:\t%f\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\n", curr_indiv_cb->index,
            										   curr_indiv_cb->patch_land_thresh,
													   root_data->current_time, 
			                                           curr_indiv_cb->current_x_pos, 
												       curr_indiv_cb->current_y_pos,
												       curr_indiv_cb->energy,
												       curr_indiv_cb->state,
													   curr_indiv_cb->ind_patch_count,
													   curr_indiv_cb->energy_gain_count,
													   curr_indiv_cb->energy_loss_count,
													   curr_indiv_cb->avg_patch_quality,
													   curr_indiv_cb->avg_search_time,
													   curr_indiv_cb->avg_distance_travelled_between_patches,
													   curr_indiv_cb->avg_distance_between_patches
													   );
				fprintf(root_data->final_output_file, "%u\t Individual %i:\t%f\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\n", 
				                                       root_data->timestamp,
													   curr_indiv_cb->index,
            										   curr_indiv_cb->patch_land_thresh,
													   root_data->current_time, 
			                                           curr_indiv_cb->current_x_pos, 
												       curr_indiv_cb->current_y_pos,
												       curr_indiv_cb->energy,
												       curr_indiv_cb->state,
													   curr_indiv_cb->ind_patch_count,
													   curr_indiv_cb->energy_gain_count,
													   curr_indiv_cb->energy_loss_count,
													   curr_indiv_cb->avg_patch_quality,
													   curr_indiv_cb->avg_search_time,
													   curr_indiv_cb->avg_distance_travelled_between_patches,
													   curr_indiv_cb->avg_distance_between_patches
													   ); 
// END: Output for data collection
			}
            
            /* Get the next CB */
            prev_indiv_cb = curr_indiv_cb;
    	    curr_indiv_cb = (PSF_INDIVIDUAL *)LIST_V2_GET_NEXT(curr_indiv_cb->list_elt);        	    

            /* If previous individual's energy has gone below zero, remove from list */
            if(prev_indiv_cb->energy < 0)
            {
// START: Output for data collection
        		fprintf(root_data->final_output_file, "%u\t Individual final stats: \t%i\t%f\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\n", 
														root_data->timestamp,
														curr_indiv_cb->index,
														curr_indiv_cb->patch_land_thresh,
													    root_data->current_time, 
			                                            curr_indiv_cb->current_x_pos, 
												        curr_indiv_cb->current_y_pos,
												        curr_indiv_cb->energy,
												        curr_indiv_cb->state,
													    curr_indiv_cb->ind_patch_count,
													    curr_indiv_cb->energy_gain_count,
													    curr_indiv_cb->energy_loss_count,
													    curr_indiv_cb->avg_patch_quality,
													    curr_indiv_cb->avg_search_time,
													    curr_indiv_cb->avg_distance_travelled_between_patches,
													    curr_indiv_cb->avg_distance_between_patches
													    );
			    printf("Individual %i has died: \t%f\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\n", curr_indiv_cb->index,
														curr_indiv_cb->patch_land_thresh,
													    root_data->current_time, 
			                                            curr_indiv_cb->current_x_pos, 
												        curr_indiv_cb->current_y_pos,
												        curr_indiv_cb->energy,
												        curr_indiv_cb->state,
													    curr_indiv_cb->ind_patch_count,
													    curr_indiv_cb->energy_gain_count,
													    curr_indiv_cb->energy_loss_count,
													    curr_indiv_cb->avg_patch_quality,
													    curr_indiv_cb->avg_search_time,
													    curr_indiv_cb->avg_distance_travelled_between_patches,
													    curr_indiv_cb->avg_distance_between_patches
													    );										    
// END: Output for data collection          	           	
            	LIST_V2_DELETE_CURRENT(&root_data->individual_list_root, prev_indiv_cb->list_elt);
            	psf_indiv_term(prev_indiv_cb);
			}
		}
        
        /* Add patches as a Poisson process: use a Poisson distribution to calculate how many
		 * patches should be added in a time interval of root_data->timestep */
		/* Start by calculating a random number to compare with the cumulative distribution */
		random_no = ((double)rand())/((double)RAND_MAX);
		/* Probability of zero patches appearing in timestep */
        curr_cum_dist = exp(-root_data->timestep*root_data->patch_app_rate);
        /* Add patches until cumulative distribution is greater than randomly calculated value */
        no_patches_added = 0;
        k_factorial = 1;
        while(curr_cum_dist < random_no)
        {
        	/* Add a patch */
        	psf_patch_init(root_data);
        	root_data->patch_count++;
        	no_patches_added++;
        	k_factorial*=no_patches_added;
        	/* Add to the cumulative distirbuion the probability of no_patches_added being added */
            curr_cum_dist += pow(root_data->timestep*root_data->patch_app_rate,no_patches_added)*
			                 exp(-root_data->timestep*root_data->patch_app_rate)/k_factorial;
		}
        
        /* Terminate patches depending on patch longevity and current patch quality */
        curr_patch_cb = (PSF_PATCH *)LIST_V2_GET_FIRST(root_data->patch_list_root);
        while(curr_patch_cb != NULL)
        {
        	/* Increment patch age */
        	curr_patch_cb->age += root_data->timestep;
// START: Output for data collection
            fprintf(root_data->patch_output_file, "%u\tPatch:\t%u\t%f\t%f\t%f\t%f\t%f\t%f\n", 
												   root_data->timestamp,
												   curr_patch_cb->index,
												   root_data->current_time, 
			                                       curr_patch_cb->x_pos,
                                                   curr_patch_cb->y_pos,
                                                   curr_patch_cb->quality,
			                                       curr_patch_cb->age,
			                                       curr_patch_cb->longevity
												   );
// END: Output for data collection
        	
            /* Get the next CB */
            prev_patch_cb = curr_patch_cb;
    	    curr_patch_cb = (PSF_PATCH *)LIST_V2_GET_NEXT(curr_patch_cb->list_elt);        	    

            /* If previous patch's energy is below zero or its age is greater than its longevity 
			 * then remove from list */
            if((prev_patch_cb->quality < 0) || (prev_patch_cb->age > prev_patch_cb->longevity))
            {
            	LIST_V2_DELETE_CURRENT(&root_data->patch_list_root, prev_patch_cb->list_elt);
            	psf_patch_term(prev_patch_cb);
			}
        }
    }
    
    /* Return */
    return;
}

/***************************************************************************************
 * Name: psf_move_individual()
 *
 * Purpose: Move the individual.
 *
 * Parameters: IN/OUT individual - the individual
 *             IN     root_data - pointer to a CB containing all the root data for PSF.
 *
 * Returns: Nothing.
 *
 * Operation: If searching:
 *             - draw a new location based on step length and cognitive map
 *             - if the line between the current and new location passes within the
 *               sensing radius of a patch that the bird is not in, and the patch is 
 *               deemed sufficiently high quality, then
 *                - move to the patch
 *                - update the state to foraging
 *                - update the patch associated to this individual
 *               else 
 *                - move to the new location
 *             - deprecate energy lost due to flying
 *            If foraging:
 *             - check the patch is still there and if not, switch to searching; else:
 *               - increase the energy of the individual due to foraging
 *               - deprecate the energy of the patch             
 *               - decide whether the individual should switch to searching, based on 
 *                 whether the quality of the patch has gone below the point where 
 *                 the time to leave is optimal (as determined by a single individual
 *                 alone seeking to minimise the probability of losing energy over time:
 *                 see Jeffries et al.).
 *            In either case, deprecate energy due to metabolic functioning
 ***************************************************************************************/
void psf_move_individual(PSF_INDIVIDUAL *individual, PSF_ROOT_DATA *root_data)
{
	/* Local variables */
	PSF_PATCH *curr_patch;
	PSF_PATCH *found_patch;
	unsigned short sample_no;
	double random_angle;
	double random_sl;
	double random_no;
	double trial_x[PSF_NO_SAMPLES];
	double trial_y[PSF_NO_SAMPLES];
	double new_x;
	double new_y;
	bool new_x_assigned;
	bool new_y_assigned;
	double weights[PSF_NO_SAMPLES];
	double sum_weights;
	bool cdm_on;
	double target_angle;
	unsigned long curr_x;
	unsigned long curr_y;
	double search_time;
	double prev_avg_search_time;
	double prev_avg_patch_quality;
	unsigned long cog_x;
	unsigned long cog_y;
//	double opt_foraging_time;
	double lambert_w_contents;
//	double foraging_start_time;
	double foraging_time;
//	unsigned long prev_patch_x;
//	unsigned long prev_patch_y;
//	double prev_patch_quality;
    double distance_to_prev_patch;
    double prev_avg_distance_between_patches;
    double prev_avg_distance_travelled_between_patches;
    double prev_x;
	double prev_y;
//	unsigned short cog_map_loc;
//	double sum_cog_map;
//	bool target_x_assigned;
//	bool target_y_assigned;
    double target_loc;

    /* Sanity checks.  If individual is searching then the patch should be NULL and time 
	 * left in patch should be zero. */
	assert(root_data != NULL);
	assert(individual != NULL);
// START: Output for testing
    //printf("individual->state:\t%i\n", individual->state);
    //printf("individual->patch:\t%i\n", individual->patch);	 
// END: Output for testing	
	assert((individual->state == PSF_FORAGING) || (individual->patch == NULL));
	assert((individual->state == PSF_FORAGING) || (individual->state == PSF_SEARCHING));
	
// START: Output for testing
//    printf("individual->state:\t%i\n", individual->state); 
// END: Output for testing
	if(individual->state == PSF_SEARCHING)
	{
		/* Save previous location for distance travelled calculation */
		prev_x = individual->current_x_pos;
		prev_y = individual->current_y_pos;
// START: Output for testing
//    printf("prev_x:\t%f\n", prev_x);
//    printf("prev_y:\t%f\n", prev_y);	 
// END: Output for testing
		
		/* Only allow cognitive direction multiplier to influence step selection function
		* after a minimum number of patches visited to avoid total cognitive map = 0 and 
		* target location remaining fixed at the 1 known location */
		cdm_on = 0;
		if(individual->ind_patch_count > 10)
		{
			cdm_on = 1;
		}
		/* Calculate target bearing between current location and target location */
		target_angle = atan2(individual->target_y-individual->current_y_pos, individual->target_x-individual->current_x_pos);	
		
		/* Find PSF_NO_SAMPLES random bearings between 0 and 2*pi and step lengths drawn from
		 * an exponential distribution with mean root_data->mean_step_length */
		sum_weights = 0;
		for(sample_no = 0; sample_no < PSF_NO_SAMPLES; sample_no++)
		{
			trial_x[sample_no] = -1;
			trial_y[sample_no] = -1;
			
			/* Cannot have a location outside the box */
			while((trial_x[sample_no] < 0) || (trial_y[sample_no] < 0) ||
			      (trial_x[sample_no] > root_data->box_width) || (trial_y[sample_no] > root_data->box_height))
			{
    			/* Number picked uniformly at random from [0,2*pi) */
	    	    random_angle = ((double)rand())*2*PSF_PI/((double)RAND_MAX + 1); // Changed so 2pi no longer repeated. Doesn't need adjusting.
		        /* This is the inverse of the cumulative exponential distribution with rate parameter 1/mean_step_length */
		        random_sl = (0-log(1-((double)rand())/((double)RAND_MAX + 1)))*root_data->mean_step_length;
		        trial_x[sample_no] = individual->current_x_pos+random_sl*cos(random_angle);
		        trial_y[sample_no] = individual->current_y_pos+random_sl*sin(random_angle);
			}
			
		    /* The numerator of a movement kernel proportional to exp(-lambda*|x0-x|+alpha*C(x)+beta*cos(theta_x-theta_x0)) where 
			 * |x0-x| is the sample step length, lambda=1/mean_step_length, C(x) is the cognitive mapat x, alpha is cog_map_mult,   
			 * theta_x0 is the target bearing, theta_x is the sample bearing and beta is cog_dir_mult */
		    weights[sample_no] = exp(root_data->cog_map_mult*individual->cog_map[(unsigned long)floor(trial_x[sample_no]) +
									    (unsigned long)(root_data->box_width)*(unsigned long)floor(trial_y[sample_no])] +
									    cdm_on*root_data->cog_dir_mult*cos(target_angle-random_angle));
			sum_weights += weights[sample_no];
// START: Output for testing
            if(isfinite(sum_weights) == 0)
            {
            	printf("Time:\t%f\tsample_no: \t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
			                    root_data->current_time,
								sample_no, 
								trial_x[sample_no], 
								trial_y[sample_no], 
								random_angle, 
								random_sl, 
								weights[sample_no], 
								sum_weights, 
								individual->cog_map[(unsigned long)floor(trial_x[sample_no])+(unsigned long)(root_data->box_width)*(unsigned long)floor(trial_y[sample_no])]); 
			}
// END: Output for testing			
// START: Output for testing
            /*printf("sample_no: \t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
			                    sample_no, 
								trial_x[sample_no], 
								trial_y[sample_no], 
								random_angle, 
								random_sl, 
								weights[sample_no], 
								sum_weights, 
								individual->cog_map[(unsigned long)floor(trial_x[sample_no])+(unsigned long)(root_data->box_width)*(unsigned long)floor(trial_y[sample_no])]); */
// END: Output for testing
		}
		
		/* Pick a random number uniformly between 0 and the sum of the weights of moving
		 * there, given the step length and the turning angle */
		random_no = ((double)rand())*sum_weights/((double)RAND_MAX); // This random_no is 
// START: Output for testing
        //printf("random_no:\t%f\n", random_no); 
    	//printf("sum_weights:\t%f\n", sum_weights); 
// END: Output for testing
		 
		/* If this random number is between the sum of the first i weights and the i+1-th weight
		 * then new position is the position given by the i-th sample */
		sum_weights = 0;
		new_x_assigned = false;
		new_y_assigned = false;
		for(sample_no = 0; sample_no < PSF_NO_SAMPLES; sample_no++)
		{
			sum_weights += weights[sample_no];
			if(sum_weights >= random_no)
			{
				new_x = trial_x[sample_no];
				new_x_assigned = true;
				new_y = trial_y[sample_no];
				new_y_assigned = true;
				break;
			}									
		}
		if(new_x_assigned == false)
		{
			fprintf(stderr,"New location cannot be chosen\n");
// START: Output for testing
        	printf("random_no:\t%f\n", random_no); 
    		printf("sum_weights:\t%f\n", sum_weights); 
// END: Output for testing
		}
        /* Sanity checks.  A sample should have been chosen. */
		assert(new_x_assigned == true);
		assert(new_y_assigned == true);
// START: Output for testing
		/*if(sum_weights <= random_no)
		{
			printf("Time:\t%f\tsum_weights\t%f\t is less than random number\t%f\n", 
									root_data->current_time, sum_weights, random_no); 
		}*/				
// END: Output for testing

// START: Output for testing
        //printf("new(x,y):\t%f\t%f\n", new_x,new_y); 
// END: Output for testing

        /* If individual has passed withing the sensing radius of the target location, choose
         * a new target location (addition as of 07/12/22, suffix of tlv - target location
		 * variable as opposed to tlf - target location fixed) */
        if((psf_dist_seg_to_point(individual->target_x,
				                    individual->target_y,
								    individual->current_x_pos,
									individual->current_y_pos,
									new_x,
									new_y) < root_data->sensing_radius))
		{
			psf_new_target_loc(individual, root_data);
// START: Output for testing
//            printf("Achieved sr of target at\t%f\t%f\tNew target\t%f\t%f\n", new_x, new_y, individual->target_x, individual->target_y); 
// END: Output for testing
		}		 
		 
		/* Loop through patches, finding out whether the line between the current position and 
		 * the proposed new one is within root_data->sensing_radius of the patch, and picking
		 * the highest quality patch if more than one are found */
		found_patch = NULL;
        if(!LIST_V2_EMPTY(root_data->patch_list_root))
        {
            curr_patch = (PSF_PATCH *)LIST_V2_GET_FIRST(root_data->patch_list_root);
            while(curr_patch != NULL)
            {
// START: Output for testing
                /*printf("Found patch at\t%i\n", curr_patch); 
                printf("Distance to\t%f\t%f\n", curr_patch->x_pos, curr_patch->y_pos);
                printf("from line between\t%f\t%f\n", individual->current_x_pos, individual->current_y_pos);
                printf("and\t%f\t%f\n", new_x, new_y);
                printf("is\t%f\n",psf_dist_seg_to_point(curr_patch->x_pos,
				                           curr_patch->y_pos,
										   individual->current_x_pos,
										   individual->current_y_pos,
										   new_x,
										   new_y)); */
// END: Output for testing
            	if((psf_dist_seg_to_point(curr_patch->x_pos,
				                          curr_patch->y_pos,
								    	  individual->current_x_pos,
										  individual->current_y_pos,
										  new_x,
										  new_y) < root_data->sensing_radius) &&
            	   ((found_patch == NULL) || (found_patch->quality < curr_patch->quality)) &&
				   ((curr_patch->x_pos != individual->prev_patch_x) || (curr_patch->y_pos != individual->prev_patch_y)))
            	{
            		/* This patch close enough to sense and better than any found previously: 
					 * set it to found_patch */
// START: Output for testing
/*                    printf("Set found_patch to\t%i at\t%f\t%f\n", curr_patch, curr_patch->x_pos, curr_patch->y_pos);
                    printf("found_patch:\t%f\t%f\t%f\t%f\t%f\n", curr_patch->x_pos,
                                                   curr_patch->y_pos,
			                                       curr_patch->age,
			                                       curr_patch->longevity,
										           curr_patch->quality); */
// END: Output for testing
            		found_patch = curr_patch;	
				}
    	        curr_patch = (PSF_PATCH *)LIST_V2_GET_NEXT(curr_patch->list_elt);
            }
	    }
	    
		if((found_patch != NULL) && (found_patch->quality > individual->patch_land_thresh))
		{
    	    /* A patch has been found that is of sufficient quality to land in.
	    	 * Move individual to this patch, associate the patch to the individual and switch 
		     * state to FORAGING */
// START: Output for testing
//            printf("Move individual to found_patch\t%i\t at\t%f\t%f\n",found_patch,found_patch->x_pos,found_patch->y_pos); 
// END: Output for testing
		    individual->current_x_pos = found_patch->x_pos;
		    individual->current_y_pos = found_patch->y_pos;
			individual->patch = found_patch;
			individual->state = PSF_FORAGING;
			individual->ind_patch_count += 1;
// START: Output for testing
//			printf("ind_patch_count:\t%i\n", individual->ind_patch_count);
// END: Output for testing			
			individual->foraging_start_time = root_data->current_time;
// START: Output for testing
//			printf("foraging_start_time:\t%f\n", individual->foraging_start_time);
// END: Output for testing
            
            /* Decide if individual gained or lost energy in prev forage-search event */
			if(individual->ind_patch_count > 1) 
			{
				if(individual->patch_start_energy > individual->energy - root_data->flying_energy_loss*root_data->timestep - root_data->timestep)
				{
					individual->energy_loss_count++;
// START: Output for testing
//			        printf("energy lost in f-s event\t%f\t%f\n", individual->patch_start_energy, individual->energy - root_data->flying_energy_loss*root_data->timestep - root_data->timestep);
//			        printf("individual->energy_loss_count\t%i\n", individual->energy_loss_count);
// END: Output for testing
				}
				else
				{
					individual->energy_gain_count++;
// START: Output for testing
//			        printf("energy gained in f-s event\t%f\t%f\n", individual->patch_start_energy, individual->energy - root_data->flying_energy_loss*root_data->timestep - root_data->timestep);
//			        printf("individual->energy_gain_count\t%i\n", individual->energy_gain_count);
// END: Output for testing
				}
			}
			individual->patch_start_energy = individual->energy - root_data->flying_energy_loss*root_data->timestep - root_data->timestep;
			// note this includes energy that will be lost after this step due to code setup.	
// START: Output for data collection
            fprintf(root_data->fs_energy_file, "Finish searching:\t%f\n", individual->patch_start_energy);
// END: Output for data collection

            /* Update average search time between patches */		
			search_time = root_data->current_time - individual->search_start_time;
// START: Output for testing
//			printf("This search\t%f\t%f\t%f\t%i\n", individual->search_start_time, root_data->current_time, search_time, individual->ind_patch_count);
// END: Output for testing
			prev_avg_search_time = individual->avg_search_time;
// START: Output for testing
//			printf("prev_avg_search_time:\t%f\n", prev_avg_search_time);
// END: Output for testing
			individual->avg_search_time = prev_avg_search_time + (search_time - prev_avg_search_time) / individual->ind_patch_count;
// START: Output for testing
//			printf("New average search time:\t%f\n", individual->avg_search_time);
// END: Output for testing

            /* Update average patch quality */
			prev_avg_patch_quality = individual->avg_patch_quality;
// START: Output for testing
//			printf("prev_avg_patch_quality:\t%f\n", prev_avg_patch_quality);
// END: Output for testing				
			individual->avg_patch_quality = prev_avg_patch_quality + (found_patch->quality - prev_avg_patch_quality) / individual->ind_patch_count;
// START: Output for testing
//			printf("New average patch quality:\t%f\n", individual->avg_patch_quality);
// END: Output for testing

            /* Update average straight line distance between patches */
            /* If this is the first patch, set prev locs to initial position */
            if(individual->ind_patch_count == 1)
            {
            	individual->prev_patch_x = individual->initial_x;
            	individual->prev_patch_y = individual->initial_y;           	
			}
            distance_to_prev_patch = sqrt(pow(individual->prev_patch_x-individual->current_x_pos,2) + pow(individual->prev_patch_y-individual->current_y_pos,2));
// START: Output for testing
/*			printf("distance_to_prev_patch\t%f\t%f\t%f\t%f\t%f\n", distance_to_prev_patch, 
			                                                    individual->prev_patch_x, 
																individual->current_x_pos, 
																individual->prev_patch_y, 
																individual->current_y_pos); */
// END: Output for testing
            prev_avg_distance_between_patches = individual->avg_distance_between_patches;
// START: Output for testing
//			printf("prev_avg_distance_between_patches:\t%f\n", prev_avg_distance_between_patches);
// END: Output for testing	
			individual->avg_distance_between_patches = prev_avg_distance_between_patches 
			                                            + (distance_to_prev_patch - prev_avg_distance_between_patches) / individual->ind_patch_count;
// START: Output for testing
//			printf("New average distance between patches:\t%f\n", individual->avg_distance_between_patches);
// END: Output for testing

		    /* Update distance travelled since patch */		
// START: Output for testing
//   		printf("distance prev travelled since patch:\t%f\n", individual->distance_travelled_since_patch);
//		    printf("old/new locations (x,y):\t%f\t%f\t%f\t%f\n", prev_x, prev_y, individual->current_x_pos, individual->current_y_pos);
// END: Output for testing  		
		    individual->distance_travelled_since_patch +=  sqrt(pow(individual->current_x_pos-prev_x,2) + pow(individual->current_y_pos-prev_y,2));
// START: Output for testing
//		    printf("distance travelled since patch:\t%f\n", individual->distance_travelled_since_patch);
// END: Output for testing

            /* Update average distance travelled between patches 
			Note: this is also in the code later when no patch is found */
			prev_avg_distance_travelled_between_patches = individual->avg_distance_travelled_between_patches;
// START: Output for testing
//			printf("prev_avg_distance_travelled_between_patches:\t%f\n", prev_avg_distance_travelled_between_patches);
// END: Output for testing	        
			individual->avg_distance_travelled_between_patches = prev_avg_distance_travelled_between_patches 
			                                            + (individual->distance_travelled_since_patch - prev_avg_distance_travelled_between_patches) 
														/ individual->ind_patch_count;	
// START: Output for testing
//			printf("New average distance travelled between patches:\t%f\n", individual->avg_distance_travelled_between_patches);
// END: Output for testing													
			individual->distance_travelled_since_patch = 0;					
// START: Output for testing
//		    printf("distance travelled since patch:\t%f\n", individual->distance_travelled_since_patch);
// END: Output for testing
			
			/* Calculate optimal time for individual to leave found_patch to optimise
			 * rate of energy gain according to Jeffries et al. model. Also corresponding
			 * patch quality leaving threshold. */
			lambert_w_contents = ((root_data->flying_energy_loss*individual->avg_search_time/found_patch->quality)-1)
									*exp(-individual->avg_search_time-1);
// START: Output for testing
//			printf("beta, rho0, 1/Lambda:\t%f\t%f\t%f\n", root_data->flying_energy_loss, found_patch->quality, individual->avg_search_time);
//			printf("lambert_w_contents:\t%.50f\n", lambert_w_contents);
// END: Output for testing

			if(lambert_w_contents>=0)
			{
				/* If Lambert W contents are greater than 0, then beta*1/Lambda > rho0
				 * i.e. the energy cost of searching is higher than the available forage in the patch. 
				 * Jeffries et al. model says maximum rate does not exist. Rate increases to -1 as t->+inf,
				 * i.e. the individual should not leave the patch. 
				 * For this model we set the opt_foraging_time as patch longevity, as this is an upper bound 
				 * for the lifespan of the patch, and the individual will not leave until it is forced to. */
				individual->opt_foraging_time = found_patch->longevity; 
				individual->curr_patch_leaving_thres = found_patch->quality*exp(-found_patch->longevity);
			}
			else
			{
				individual->opt_foraging_time = -psf_lambert_w(lambert_w_contents, root_data)-individual->avg_search_time-1;
// START: Output for testing				
//				printf("psf_lambert_w(lambert_w_contents):\t%f\n", psf_lambert_w(lambert_w_contents, root_data));	
// END: Output for testing
				individual->curr_patch_leaving_thres = found_patch->quality*exp(-individual->opt_foraging_time);				
			}			
// START: Output for testing
//			printf("opt_foraging_time:\t%f\n", individual->opt_foraging_time);
//			printf("curr_patch_leaving_thres:\t%f\n", individual->curr_patch_leaving_thres);
// END: Output for testing	
			
// START: Output for testing
//			printf("prev patch details for cog map:\t%i\t%i\t%f\n", individual->prev_patch_x, individual->prev_patch_y, individual->prev_patch_quality);
// END: Output for testing				
     		/* Update cognitive map: increment the cognitive map for every grid point within the sensing 
			 * radius of the previous patch by the quality of the patch */
      		for(curr_y = max(0,(unsigned long)floor(individual->prev_patch_y)-root_data->sensing_radius); 
			    curr_y < min(root_data->box_height,(unsigned long)floor(individual->prev_patch_y)+root_data->sensing_radius+1);
				curr_y++)
			{
          		for(curr_x = max(0,(unsigned long)floor(individual->prev_patch_x)-root_data->sensing_radius); 
			    curr_x < min(root_data->box_width,(unsigned long)floor(individual->prev_patch_x)+root_data->sensing_radius+1);
  				curr_x++)
	    		{
	    			if(sqrt(pow((float)curr_x-individual->prev_patch_x,2)+pow((float)curr_y-individual->prev_patch_y,2))<root_data->sensing_radius)
	    			{
// START: Output for testing
//                      printf("Increment cognitive map at\t%i\t%i\n",curr_x, curr_y);
//						printf("Cognitive map value at\t%i\t%i\t was\t%f\n",curr_x, curr_y, individual->cog_map[curr_x+(unsigned long)(root_data->box_width)*curr_y]); 
// END: Output for testing
	    				individual->cog_map[curr_x+(unsigned long)(root_data->box_width)*curr_y] += individual->prev_patch_quality;
// START: Output for testing
//               	    printf("Cognitive map value at\t%i\t%i\t is now\t%f\n",curr_x, curr_y, individual->cog_map[curr_x+(unsigned long)(root_data->box_width)*curr_y]); 
// END: Output for testing
					}
		    	} 			
			}
			/* Store details of current patch to update cognitive map when next patch is found */ 
			individual->prev_patch_x = found_patch->x_pos;
			individual->prev_patch_y = found_patch->y_pos;
			individual->prev_patch_quality = found_patch->quality;
// START: Output for testing
//			printf("save patch details for cog map:\t%f\t%f\t%f\n", individual->prev_patch_x, individual->prev_patch_y, individual->prev_patch_quality);
// END: Output for testing   		
		} 
        else		 
		{
    		/* No patch has been found worth landing in.  Move the individual to the newly-found location */
		    individual->current_x_pos = new_x;
		    individual->current_y_pos = new_y;
		    
		    /* Update distance travelled since patch.
			 * Note: this is also in the code prev when a patch is found */		
// START: Output for testing
//   		printf("distance prev travelled since patch:\t%f\n", individual->distance_travelled_since_patch);
//		    printf("old/new locations (x,y):\t%f\t%f\t%f\t%f\n", prev_x, prev_y, individual->current_x_pos, individual->current_y_pos);
// END: Output for testing  		
		    individual->distance_travelled_since_patch +=  sqrt(pow(individual->current_x_pos-prev_x,2) + pow(individual->current_y_pos-prev_y,2));
// START: Output for testing
//		    printf("distance travelled since patch:\t%f\n", individual->distance_travelled_since_patch);
// END: Output for testing 
		}
		  
		/* Deprecate energy due to flying */
		individual->energy -= root_data->flying_energy_loss*root_data->timestep;
		
		/* Cognitive map decay */
		for(cog_y = 0; cog_y < root_data->box_height; cog_y++)
		{
			for(cog_x = 0; cog_x < root_data->box_height; cog_x++)
			{
// START: Output for testing
				//printf("Cognitive map value at\t%i\t%i\t was\t%f\n",cog_x, cog_y, individual->cog_map[cog_x+(unsigned long)(root_data->box_width)*cog_y]); 
// END: Output for testing
				individual->cog_map[cog_x+(unsigned long)(root_data->box_width)*cog_y] -= 
				                                                   individual->cog_map[cog_x+(unsigned long)(root_data->box_width)*cog_y]
																   *root_data->timestep
																   *(root_data->cog_map_decay);
// START: Output for testing
                //printf("Cognitive map value at\t%i\t%i\t is now\t%f\n",cog_x, cog_y, individual->cog_map[cog_x+(unsigned long)(root_data->box_width)*cog_y]); 
// END: Output for testing
			}
		}	
	}
	else
	{
		/* Individual is foraging */
 		assert(individual->state == PSF_FORAGING);
 		if(individual->patch == NULL)
 		{
     		/* The patch is no longer there.  Switch to searching */
     		individual->state = PSF_SEARCHING;
     		
     		/* Store time individual leaves patch */
    		individual->search_start_time = root_data->current_time;
// START: Output for testing
//            printf("Individual->search_start_time\t%f\n",individual->search_start_time);
// END: Output for testing
		}
 		else
 		{
     		assert(individual->patch != NULL);
     		/* Increase energy of individual due to foraging according to Jeffries et al. model */
     		individual->energy += (individual->patch)->quality*root_data->timestep;
     		
     		/* Deprecate energy of patch, according to Jeffries et al. model */
     		(individual->patch)->quality -= (individual->patch)->quality*root_data->timestep;
     		
     		/* Calculate current time spent in patch */
     		foraging_time = root_data->current_time-individual->foraging_start_time;
// START: Output for testing
//          printf("foraging_time:\t%f\n",foraging_time);
// END: Output for testing
     		
     		if((individual->curr_patch_leaving_thres > (individual->patch)->quality) || ((individual->patch)->age >= (individual->patch)->longevity))
     		{
// START: Output for testing
//			    printf("foraging_start_time:\t%f\n",individual->foraging_start_time);
//              printf("individual->foraging_time:\t%f\n",foraging_time);
//              printf("individual->opt_foraging_time:\t%f\n",individual->opt_foraging_time);
//				printf("(Individual->patch)->age\t%f\n",(individual->patch)->age);
// 	          	printf("(Individual->patch)->longevity\t%f\n",(individual->patch)->longevity);
// END: Output for testing
    			individual->state = PSF_SEARCHING;
    			individual->patch = NULL; 

    			/* Store time individual leaves patch */
    			individual->search_start_time = root_data->current_time;
// START: Output for testing
//        		printf("Individual->search_start_time\t%f\n",individual->search_start_time);
// END: Output for testing
// START: Output for data collection
            	fprintf(root_data->fs_energy_file, "Finish foraging:\t%f\n", individual->energy - root_data->timestep);
            	// note this includes energy that will be lost after this step due to code setup.
// END: Output for data collection
				
				/* Choosing next target location target_x, target_y */
				psf_new_target_loc(individual, root_data);			
				 		
// Old code kept for reference purposes			     		
     		/* When time spent in patch exceeds optimal foraging time tau_f^* 
			* individual switches to searching, as in Jeffries et. al model 
     		if(foraging_time >= individual->opt_foraging_time)
     		{
// START: Output for testing
			    printf("foraging_start_time:\t%f\n",individual->foraging_start_time);
                printf("individual->foraging_time:\t%f\n",foraging_time);
                printf("individual->opt_foraging_time:\t%f\n",individual->opt_foraging_time);
// END: Output for testing
    			individual->state = PSF_SEARCHING;
    			individual->patch = NULL; */
    			
     		
     		/* If patch quality goes below 1, switch to searching.  This is optimal for a single lone
			 * individual seeking to minimise the risk of unsuccessful foraging, in dimensionless units 
     		if((individual->patch)->quality < 1)
     		{
    			individual->state = PSF_SEARCHING;
    			individual->patch = NULL; */
    										
			}
		}
	}
	/* Deprecate energy due to metabolic needs (this is delta-t times 1, in dimensionless units) */
	individual->energy -= root_data->timestep;
	
    /* Return */
    return;
}

/***************************************************************************************
 * Name: psf_dist_seg_to_point()
 *
 * Purpose: Find the distance from a line segment to a point.
 *
 * Parameters: IN     point_x - x-value of the point 
 *             IN     point_y - y-value of the point 
 *             IN     line_pt1_x - x-value of one end of the line segment
 *             IN     line_pt1_y - y-value of one end of the line segment
 *             IN     line_pt2_x - x-value of the other end of the line segment
 *             IN     line_pt2_y - y-value of the other end of the line segment
 *
 * Returns: distance - Distance from the point to the line segment
 *
 * Operation: Find the shortest distance from the (point_x, point_y) to the line extended 
 *            indefinitely from either end.  Fine the point (a,b) on the extended line where 
 *            the point (point_x,point_y) is closest.  If (a,b) is on the line segment 
 *            then return the distance from (a,b) to (point_x, point_y).  If not, return 
 *            the shortest distance between the distance from (point_x, point_y) to 
 *            (line_pt1_x,line_pt1_y) and the distance from (point_x, point_y) to 
 *            (line_pt2_x,line_pt2_y).
 ***************************************************************************************/
double psf_dist_seg_to_point(double point_x, 
                             double point_y, 
							 double line_pt1_x, 
							 double line_pt1_y,
							 double line_pt2_x, 
							 double line_pt2_y)
{
	/* Local variables */
	double distance;
	double dist_whole_line;
	double x_perp_drop;
	double y_perp_drop;
	
    /* Sanity checks. */
    
    /* Find the point where the perpendicular from (point_x, point_y) meets the extended line */
	/* x_perp_drop = ((line_pt1_y-point_y)*(line_pt2_x-line_pt1_x)*(line_pt2_y-line_pt1_y)+
	               line_pt1_x*pow(line_pt2_y-line_pt1_y,2)+point_x*pow(line_pt2_x-line_pt1_x,2))/
	               (pow(line_pt2_y-line_pt1_y,2)+pow(line_pt2_x-line_pt1_x,2));
	y_perp_drop = ((line_pt1_x-point_x)*(line_pt2_y-line_pt1_y)*(line_pt2_x-line_pt1_x)+
	               line_pt1_y*pow(line_pt2_x-line_pt1_x,2)+point_y*pow(line_pt2_y-line_pt1_y,2))/
	               (pow(line_pt2_x-line_pt1_x,2)+pow(line_pt2_y-line_pt1_y,2));  */
	x_perp_drop = ((point_y-line_pt1_y)*(line_pt2_x-line_pt1_x)*(line_pt2_y-line_pt1_y)+
	               line_pt1_x*pow(line_pt2_y-line_pt1_y,2)+point_x*pow(line_pt2_x-line_pt1_x,2))/
	               (pow(line_pt2_y-line_pt1_y,2)+pow(line_pt2_x-line_pt1_x,2));
	y_perp_drop = ((point_x-line_pt1_x)*(line_pt2_y-line_pt1_y)*(line_pt2_x-line_pt1_x)+
	               line_pt1_y*pow(line_pt2_x-line_pt1_x,2)+point_y*pow(line_pt2_y-line_pt1_y,2))/
	               (pow(line_pt2_x-line_pt1_x,2)+pow(line_pt2_y-line_pt1_y,2)); 
// START: Output for testing
//    printf("x_perp_drop, y_perp_drop \t%f\t%f\n", x_perp_drop, y_perp_drop); 
// END: Output for testing
	
	/* Check if (x_perp_drop, y_perp_drop) is between (line_pt1_x,line_pt1_y) and (line_pt2_x,line_pt2_y) */
	if(((line_pt2_x-x_perp_drop)*(x_perp_drop-line_pt1_x) >= 0) && 
	   ((line_pt2_y-y_perp_drop)*(y_perp_drop-line_pt1_y) >= 0))
	{
		/* (x_perp_drop, y_perp_drop) is between (line_pt1_x,line_pt1_y) and (line_pt2_x,line_pt2_y) so return
		 * the distance from (point_x, point_y) to the extended line */
        distance = PSF_DIST_LINE_TO_POINT(point_x,point_y,line_pt1_x,line_pt1_y,line_pt2_x,line_pt2_y);
        
		/*distance = fabs((point_x*(line_pt2_y-line_pt1_y)-point_y*(line_pt2_x-line_pt1_x)
		*            +line_pt2_x*line_pt1_y-line_pt1_x*line_pt2_y))/
		*			(sqrt(pow(line_pt2_x-line_pt1_x,2)+pow(line_pt2_y-line_pt1_y,2))); */
// START: Output for testing
        //printf("distance \t%f\n", distance); 
        //printf("perp point check \t%f\t%f\n", (line_pt2_x-x_perp_drop)*(x_perp_drop-line_pt1_x), (line_pt2_y-y_perp_drop)*(y_perp_drop-line_pt1_y)); 
        //printf("inputs psf_dist \t%f\t%f\t%f\t%f\t%f\t%f\n", point_x,point_y,line_pt1_x,line_pt1_y,line_pt2_x,line_pt2_y);
// END: Output for testing
	}
	else
	{
		/* (x_perp_drop, y_perp_drop) is not between (line_pt1_x,line_pt1_y) and (line_pt2_x,line_pt2_y) so return
		 * the distance from (point_x, point_y) to either (line_pt1_x,line_pt1_y) or (line_pt2_x,line_pt2_y), whichever
		 * is lower */
		distance = min(sqrt(pow(point_x-line_pt1_x,2)+pow(point_y-line_pt1_y,2)),
		               sqrt(pow(point_x-line_pt2_x,2)+pow(point_y-line_pt2_y,2))); 
	}
	
    /* Return */
    return(distance);
}

/***************************************************************************************
 * Name: psf_lambert_w()
 *
 * Purpose: Calculate the output value of the Lambert W function.
 *
 * Parameters: IN     lambert_w_input 
 *
 * Returns: lambert_w_output
 *
 * Operation: 
 *
 ***************************************************************************************/

double psf_lambert_w(double lambert_w_input, PSF_ROOT_DATA *root_data)

{
	/* Local variables */
	double lambert_w_output;
	unsigned short input_no;
	double sample_input;
	double sample_diff;
	double min_diff = 10;
	double min_diff_input;
	
// START: Output for testing
//	printf("smallest psf_lambert_w input:\t%.30f\n", root_data->lambert_w[2*(root_data->no_lambert_w_inputs-1)]);
// END: Output for testing
	if((root_data->lambert_w[2*(root_data->no_lambert_w_inputs-1)] < lambert_w_input))
	{
		printf("%u\tTime:\t%f\tWarning: psf_lambert_w input too small\t%.60f\n",root_data->timestamp, root_data->current_time, lambert_w_input);
	}
		
	
	
	for(input_no = 0; input_no < root_data->no_lambert_w_inputs; input_no++)
	{
		sample_input = root_data->lambert_w[2*input_no];
		sample_diff = fabs(sample_input-lambert_w_input);
// START: Output for testing
		//printf("psf_lambert_w sampling:\t%i\t%f\t%.30f\t%f\t%f\n", input_no, lambert_w_input, sample_input, sample_diff, min_diff);
// END: Output for testing
		if((min_diff > sample_diff)) 
		{
			min_diff = sample_diff;
			min_diff_input = sample_input;
			lambert_w_output = root_data->lambert_w[2*input_no+1];
// START: Output for testing
			//printf("psf_lambert_w new min_diff::\t%i\t%f\t%.30f\t%f\n", input_no, min_diff, min_diff_input, lambert_w_output);
// END: Output for testing
		}
	}
// START: Output for testing
//	printf("psf_lambert_w:\t%f\t%f\t%f\n", lambert_w_input, min_diff_input, lambert_w_output);
// END: Output for testing	

	/* Return */	
	return(lambert_w_output);
}

/***************************************************************************************
 * Name: psf_new_target_loc()
 *
 * Purpose: Choose a new target location based on individual's cognitive map
 *
 * Parameters: IN     lambert_w_input 
 *
 * Returns: target_loc
 *
 * Operation: 
 *
 ***************************************************************************************/

void psf_new_target_loc(PSF_INDIVIDUAL *individual, PSF_ROOT_DATA *root_data)
/* Choosing next target location target_x, target_y */
{
    /* Local variables */
	unsigned short cog_map_loc;
	double sum_cog_map;
	double random_no;
	bool target_assigned;
								
	/* Calculate sum_cog_map */				
	sum_cog_map = 0;
	for(cog_map_loc = 0; cog_map_loc < (root_data->box_width*root_data->box_height); cog_map_loc++)
	{
		sum_cog_map += individual->cog_map[cog_map_loc];
	}
	/* Chose random number between 0 and sum_cog_map */					
	random_no = ((double)rand())*sum_cog_map/((double)RAND_MAX);
// START: Output for testing
//  printf("random_no, sum_cog_map:\t%f\t%f\n",random_no, sum_cog_map);
// END: Output for testing
	/* Assign target_x and target_y */			
	sum_cog_map = 0;
	target_assigned = false;
	for(cog_map_loc = 0; cog_map_loc < (root_data->box_width*root_data->box_height); cog_map_loc++)
	{
		sum_cog_map += individual->cog_map[cog_map_loc];
		if(sum_cog_map >= random_no)
		{
// START: Output for testing
//            printf("random_no, sum_cog_map:\t%f\t%f\n", random_no, sum_cog_map);
// END: Output for testing
            individual->target_y = floor(cog_map_loc/root_data->box_width);				        
			individual->target_x = cog_map_loc - (individual->target_y*root_data->box_width);
			target_assigned = true;
			break;
		}									
	}		        
// START: Output for testing
/*    printf("target_loc found\t%f\t%i\t%i\t%i\t%i\t%f\n", root_data->current_time, 
										(unsigned long)individual->current_x_pos , 
										(unsigned long)individual->current_y_pos, 
										(unsigned long)individual->target_x, 
										(unsigned long)individual->target_y, 
										individual->cog_map[(unsigned long)(individual->target_x + (root_data->box_width)*individual->target_y)]); */
// END: Output for testing
    /* Sanity checks.  A sample should have been chosen. */
	if(target_assigned == false)
	{
		fprintf(stderr,"New target location cannot be chosen\n");
	}
	assert(target_assigned == true);
	
	/* return */
    return;	
}


