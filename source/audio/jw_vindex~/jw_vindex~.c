/*
Copyright (c) 2022.  The Regents of the University of California (Regents).
All Rights Reserved.

Permission to use, copy, modify, and distribute this software and its
documentation for educational, research, and not-for-profit purposes, without
fee and without a signed licensing agreement, is hereby granted, provided that
the above copyright notice, this paragraph and the following two paragraphs
appear in all copies, modifications, and distributions.  Contact The Office of
Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley,
CA 94720-1620, (510) 643-7201, for commercial licensing opportunities.

Written by Jeremy L. Wagner, The Center for New Music and Audio Technologies,
University of California, Berkeley.

     IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
     SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
     ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
     REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

     REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
     LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
     FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING
     DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".
     REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
     ENHANCEMENTS, OR MODIFICATIONS.
*/
//
//  jw_vindex~.c
//  vindex_tilde
//  This is a sandbox attempt at returning a vector from a buffer at some specified (int) index.
//  Load that buffer into FFTW, perform FFT, then do some peak finding and rough estimation/refinement of
//  bin frequency. Output as a list of likely frequencies.
//
//  Created by Jeremy Wagner on 8/11/22.
//

#include "jw_vindex~.h" //probably not needed

#include "ext.h"
#include "ext_obex.h"
#include "ext_common.h" // contains CLAMP macro
#include "z_dsp.h"
#include "ext_buffer.h"
#include "fftw3.h"
#include "time.h"

//struct to contain analysis window
typedef struct _slice {
    double *mag_spec;
    double *phase_spec;
    double sum;
    double max_peak;
    int num_peaks;
    long *peaks;
} t_slice;

//struct for object
typedef struct _jw_vindex {
    t_pxobject l_obj;
    long l_vector;
    t_buffer_ref *l_buffer_reference;
    long l_chan;
    //t_buffer_ref *o_buffer_reference;
    //long o_chan;
    long sample_vector_size;    //length of vector we will pull from the buffer
    long cursor;                //cursor in buffer (the analysis point)
    long cursor2;               //cursor2 in buffer (the second analysis point)
    long *analysis_points;      //an array of analysis points for determining decay rates
    double *ap0;
    double *ap1;
    double *ap2;
    double *ap3;
    double *ap4;
    
    void *out;                  //outlet
    void *f_out;
    void *d_out;                //dump outlet
    long fft_size;
    fftw_plan p;
    double *in;
    double *outs;
    double *mag_spec;
    double *phase_spec;
    double thresh;
    double sum;
    double max_peak;
    int num_peaks;
    long *peaks;
    float sr;
    double *cooked;
    t_slice *slices;
    long num_samples;
    
} t_jw_vindex;



//prototypes
void jw_vindex_perform64(t_jw_vindex *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void jw_vindex_dsp64(t_jw_vindex *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void jw_vindex_int(t_jw_vindex *x, long n);
void jw_vindex_set(t_jw_vindex *x, t_symbol *s);
void jw_vindex_setvsize(t_jw_vindex *x, long n);
void jw_vindex_getvsize(t_jw_vindex *x);
void jw_vindex_set_fft_size(t_jw_vindex *x, long n);
void *jw_vindex_new(t_symbol *s, long chan);
void jw_vindex_free(t_jw_vindex *x);
t_max_err jw_vindex_notify(t_jw_vindex *x, t_symbol *s, t_symbol *msg, void *sender, void *data);
void jw_vindex_in1(t_jw_vindex *x, long n);
void jw_vindex_assist(t_jw_vindex *x, void *b, long m, long a, char *s);
void jw_vindex_dblclick(t_jw_vindex *x);
void jw_vindex_set_thresh(t_jw_vindex *x, double n);
void jw_vindex_bang(t_jw_vindex *x);
void jw_vindex_list_out(t_jw_vindex *x, double* a, long l);
void jw_vindex_list(t_jw_vindex *x, t_symbol *msg, long argc, t_atom *argv);
double *exp_fit(long *xVals, double *yVals, long n);

void hann_window(t_jw_vindex *x);

//class
static t_class *jw_vindex_class;

C74_EXPORT void ext_main(void *r)
{
    t_class *c = class_new("jw_vindex~", (method)jw_vindex_new, (method)jw_vindex_free, sizeof(t_jw_vindex), 0L, A_SYM, A_DEFLONG, 0);

    class_addmethod(c, (method)jw_vindex_dsp64, "dsp64", A_CANT, 0);
    class_addmethod(c, (method)jw_vindex_set, "set", A_SYM, 0);
    class_addmethod(c, (method)jw_vindex_in1, "in1", A_LONG, 0);
    class_addmethod(c, (method)jw_vindex_int, "int", A_LONG, 0);
    class_addmethod(c, (method)jw_vindex_assist, "assist", A_CANT, 0);
    class_addmethod(c, (method)jw_vindex_dblclick, "dblclick", A_CANT, 0);
    class_addmethod(c, (method)jw_vindex_notify, "notify", A_CANT, 0);
    class_addmethod(c, (method)jw_vindex_setvsize, "set_vector_size", A_LONG, 0);
    class_addmethod(c, (method)jw_vindex_getvsize, "get_vector_size", A_DEFSYM, 0);
    class_addmethod(c, (method)jw_vindex_set_fft_size, "fft_size", A_LONG, 0);
    class_addmethod(c, (method)jw_vindex_set_thresh, "set_thresh", A_FLOAT, 0);
    class_addmethod(c, (method)jw_vindex_bang, "bang", A_CANT, 0);                  //for experimental purposes
    class_addmethod(c, (method)jw_vindex_list, "list", A_CANT, 0);
    
    class_dspinit(c);
    class_register(CLASS_BOX, c);
    jw_vindex_class = c;
}

void jw_vindex_perform64(t_jw_vindex *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    t_double    *in = ins[0];
    t_double    *out = outs[0];
    int            n = sampleframes;
    t_float        *tab;
    double        temp;
    double        f;
    long        index, chan, frames, nc;
    t_buffer_obj    *buffer = buffer_ref_getobject(x->l_buffer_reference);

    tab = buffer_locksamples(buffer);
    if (!tab)
        goto zero;

    frames = buffer_getframecount(buffer);
    nc = buffer_getchannelcount(buffer);
    x->sr = buffer_getsamplerate(buffer);
    chan = MIN(x->l_chan, nc);
    while (n--) {
        temp = *in++;
        f = temp + 0.5;
        index = f;
        if (index < 0)
            index = 0;
        else if (index >= frames)
            index = frames - 1;
        if (nc > 1)
            index = index * nc + chan;
        *out++ = tab[index];
    }
    buffer_unlocksamples(buffer);
    return;
zero:
    while (n--)
        *out++ = 0.0;
}

static void print_result(float bw, t_jw_vindex *x) {
    for(int k=0;k<x->fft_size / 2;k++)
        post("%d: (%f) (%f, %fi), %f, %f",
             k,
             k*bw,
             x->outs[2*k],
             x->outs[2*k+1],
             sqrt(pow(x->outs[2*k],2)+pow(x->outs[2*k+1],2)),  //could remove the sqrt here
             atan2(x->outs[2*k+1],x->outs[2*k])
             );
}

//when we get an int, set the cursor and print the fft input at that point in the target buffer, then the fft of that window
//needs to be refactored
void jw_vindex_int(t_jw_vindex *x, long n)
{
    if(n>=0)
    {
        x->cursor = n;
        
        t_float *tab;
        t_buffer_obj    *buffer = buffer_ref_getobject(x->l_buffer_reference);
        x->sr = buffer_getsamplerate(buffer);
        tab = buffer_locksamples(buffer);
        if(!tab)
            goto zero;
        //get buffer length. If window at cursor exceeds buffer length, truncate window.
        long frames = buffer_getframecount(buffer);
        long i = x->cursor + x->fft_size;
        i = MIN(x->cursor, frames - x->fft_size);
        
        //get channels
        long nc = buffer_getchannelcount(buffer);
        long chan = MIN(x->l_chan, nc);
        
        //load window into fft input
        for(int j=0; j< x->fft_size;j++){
            x->in[j] = tab[(j+i)*nc+chan];
            //x->in[2*j+1] = 0;  //no imaginary component
            //post("%d: %f", j, x->in[j]);
            
        }
        
        //perform fft
        hann_window(x);
        clock_t t1, t2;
        t1=clock();
        fftw_execute(x->p);
        t2 = clock();
        post("jw_vindex: fft took %f s", (double)(t2-t1)/CLOCKS_PER_SEC);
        
        //find bin width based on window size and sample rate
        float bw = x->sr / x->fft_size;
        //print_result(bw, x);
        
        //derive magnitude & phase, find sum and max
        x->sum = 0;
        double temp=0;
        for(int i=0;i<x->fft_size;i++){
            temp = sqrt(pow(x->outs[2*i],2) + pow(x->outs[2*i+1],2));
            x->sum+=temp;
            if(temp>=x->max_peak) x->max_peak = temp;
            x->mag_spec[i]=temp;
            x->phase_spec[i]=atan2(x->outs[2*i+1],x->outs[2*i]);
        }
        
        
        
        //find peaks
        int c = 0;
        for(int i=1;i<x->fft_size-1;i++)
        {
            if((x->mag_spec[i-1]<x->mag_spec[i]) && (x->mag_spec[i+1]<x->mag_spec[i]))
            {
                //now check if the detected peak is within the threshold of our max peak
                if(20*log10(x->mag_spec[i] / x->max_peak) > x->thresh)
                {
                    //only write if it is loud enough
                    x->peaks[c]=i;
                    c++;
                }
               
            }
            x->peaks[c]=-1;
        }
        //let's track the number of peaks we've found (+1)
        x->num_peaks=c;
        c=0;
        
        
        
        
        while(x->peaks[c]>=0){
            post("jw_vindex: peak at bin %ld: (%f Hz)", x->peaks[c], x->peaks[c]*bw);
            c++;
        }
        
        //cook the pitch with a fractional bin analysis
        // /fractional_bins = [ 0, log(/spectrum[[/i+1]] / /spectrum[[/i -1]]) / (2 * log(pow(/spectrum[[/i]],2) / (/spectrum[[/i-1]] * /spectrum[[/i+1]]))), 0 ],
        for(int i=0; i<x->num_peaks;i++){
            int ind = x->peaks[i];
            double f = ind;
            //the following can be found in the literature on fractional bin extraction
            //the log functions improve accuracy, but could be omitted or other numerical methods tried for efficiency
            f += log(x->mag_spec[ind+1]/x->mag_spec[ind - 1]) / (2*log(pow(x->mag_spec[ind],2) / (x->mag_spec[ind+1] * x->mag_spec[ind - 1])) );
            x->cooked[2*i] = f*bw;      //add cooked frequency to output list
            x->cooked[2*i+1] = x->mag_spec[ind] / x->max_peak;  //add normalized amplitude to output list (for now)
            post("jw_vindex: cooked bin %f: (%f Hz)", f, f*bw);
            
        }
        jw_vindex_list_out(x, x->cooked, x->num_peaks * 2);     //list the cooked (frequency, amplitude) pairs out the outlet
        buffer_unlocksamples(buffer); //maybe move this earlier
        return;
        
        
    zero:
        outlet_float(x->f_out, 0.0);
        post("jw_vindex: Error: did not get buffer.");
    } else {
        x->cursor = 0;
        post("jw_vindex: Cursor values must be integers greater than zero. Setting to zero.");
    }
}

void jw_vindex_list(t_jw_vindex *x, t_symbol *msg, long argc, t_atom *argv){
    if(argc != 2){
        error("jw_vindex: only supports list length of 2: (cursor1, cursor2)");
        return;
    }
    
    long buffer_len = buffer_getframecount(buffer_ref_getobject(x->l_buffer_reference));
    //error("jw_vindex: buffer framecount: %ld", buffer_len);

    long c1 = atom_getlong(argv);
    long c2 = atom_getlong(argv +1);
    if(c1<0 || c1 >= buffer_len)
    {
        //buffer has no negative indices
        if(c1<0) c1=0;
        //need at least 1 fft frame of runway before end of buffer
        if(c1>buffer_len - 1 - x->fft_size) c1 = buffer_len - 1 - x->fft_size;
        post("jw_vindex: Error: cursor1 position must be between zero and buffer length. Setting to %ld", c1);
        
    } else if(c2<0 || c2 >= buffer_len)
    {
        if(c2<0) c2=0;
        if(c2>buffer_len - 1 - x->fft_size) c2 = buffer_len - 1 - x->fft_size;
        post("jw_vindex: Error: cursor2 position must be between zero and buffer length. Setting to %ld", c2);

    }
    
    x->cursor = c1;
    x->cursor2 = c2;
    post("jw_vindex: received cursor positions %ld, %ld", x->cursor, x->cursor2);  //this works
    
    //check that our cursors are in order
    if(c1>c2){
        error("jw_vindex: cursor position 2 must be greater than cursor position 1");
        return;
    }
    
    //set analysis points
    long step = (c2-c1)/4;
    x->analysis_points[0] = c1;
    x->analysis_points[1] = c1+step;
    x->analysis_points[2] = c1+2*step;
    x->analysis_points[3] = c1+3*step;
    x->analysis_points[4] = c2;
    
    
    //TODO:
    //conduct jw_vindex_int operations at c1, then perform ffts at the remaining 4 points, logging the amplitudes at bins identified as peaks at c1
    
        
        t_float *tab;
        t_buffer_obj    *buffer = buffer_ref_getobject(x->l_buffer_reference);
        x->sr = buffer_getsamplerate(buffer);
        tab = buffer_locksamples(buffer);
        if(!tab)
            goto zero;
        //get buffer length. If window at cursor exceeds buffer length, truncate window.
        long frames = buffer_getframecount(buffer);
        long i = x->cursor + x->fft_size;
        i = MIN(x->cursor, frames - x->fft_size);
        
        //get channels
        long nc = buffer_getchannelcount(buffer);
        long chan = MIN(x->l_chan, nc);
        
        
        //load window into fft input
        for(int j=0; j< x->fft_size;j++){
            x->in[j] = tab[(j+1)*nc+chan];
        }
        
        //perform fft
        hann_window(x);
        fftw_execute(x->p);

        
        //find bin width based on window size and sample rate
        float bw = x->sr / x->fft_size;
        
        //derive magnitude & phase, find sum and max
        x->sum = 0;
        double temp=0;
        for(int i=0;i<x->fft_size;i++){
            temp = sqrt(pow(x->outs[2*i],2) + pow(x->outs[2*i+1],2));
            x->sum+=temp;
            if(temp>=x->max_peak) x->max_peak = temp;
            x->mag_spec[i]=temp;
            x->phase_spec[i]=atan2(x->outs[2*i+1],x->outs[2*i]);
        }

        //find peaks
        int c = 0;
        for(int i=1;i<x->fft_size-1;i++)
        {
            if((x->mag_spec[i-1]<x->mag_spec[i]) && (x->mag_spec[i+1]<x->mag_spec[i]))
            {
                //now check if the detected peak is within the threshold of our max peak
                if(20*log10(x->mag_spec[i] / x->max_peak) > x->thresh)
                {
                    //only write if it is loud enough
                    x->peaks[c]=i;
                    c++;
                }
               
            }
            x->peaks[c]=-1;
        }
        //let's track the number of peaks we've found (+1)
        x->num_peaks=c;
        c=0;
        
        
        
        
        while(x->peaks[c]>=0){
            post("jw_vindex: peak at bin %ld: (%f Hz)", x->peaks[c], x->peaks[c]*bw);
            c++;
        }
        
        //cook the pitch with a fractional bin analysis
        // /fractional_bins = [ 0, log(/spectrum[[/i+1]] / /spectrum[[/i -1]]) / (2 * log(pow(/spectrum[[/i]],2) / (/spectrum[[/i-1]] * /spectrum[[/i+1]]))), 0 ],
        for(int i=0; i<x->num_peaks;i++){
            int ind = x->peaks[i];
            double f = ind;
            //the following can be found in the literature on fractional bin extraction
            //the log functions improve accuracy, but could be omitted or other numerical methods tried for efficiency
            f += log(x->mag_spec[ind+1]/x->mag_spec[ind - 1]) / (2*log(pow(x->mag_spec[ind],2) / (x->mag_spec[ind+1] * x->mag_spec[ind - 1])) );
            x->cooked[2*i] = f*bw;      //add cooked frequency to output list
            x->cooked[2*i+1] = x->mag_spec[ind] / x->max_peak;  //add normalized amplitude to output list (for now)
            post("jw_vindex: cooked bin %f: (%f Hz)", f, f*bw);
            
        }
        jw_vindex_list_out(x, x->cooked, x->num_peaks * 2);     //list the cooked (frequency, amplitude) pairs out the outlet
        buffer_unlocksamples(buffer); //maybe move this earlier
        return;
        
        
    zero:
        outlet_float(x->f_out, 0.0);
        post("jw_vindex: Error: did not get buffer.");
    
    
    
    //work out A,B for A*exp(-B*x) by performing exponential fitting
    //return cooked (freq, amplitude, decay_rate) triples
    
    
    
    

}




void jw_vindex_set(t_jw_vindex *x, t_symbol *s)
{
    if (!x->l_buffer_reference)
        x->l_buffer_reference = buffer_ref_new((t_object *)x, s);
    else
        buffer_ref_set(x->l_buffer_reference, s);
    
    //the buffer may have a different sample rate.  Let's find out what it is and reset our SR to match.
    t_buffer_obj    *buffer = buffer_ref_getobject(x->l_buffer_reference);
    x->sr = buffer_getsamplerate(buffer);
    
}

void jw_vindex_setvsize(t_jw_vindex *x, long n)
{
    if(n>0){
        x->sample_vector_size =n;
        post("jw_vindex: Hey, it worked! Vector size will be set to %d", n);
    } else {
        x->sample_vector_size = 1;
        post("jw_vindex: Sample Vector Size must be a positive integer. Setting to 1");
    }
}

void jw_vindex_set_fft_size(t_jw_vindex *x, long n)
{
    if(n>0){
        //bitwise-& checks for power of 2
        if((n & (n-1)) == 0){
            x->fft_size = n;
            
            //when resetting fft size we need to free and re-allocate memory
            //otherwise, we get a crash
            if(x->in != NULL) fftw_free((char *)x->in);
            if(x->outs !=NULL) fftw_free((char *)x->outs);
            x->in = (double *) fftw_malloc(sizeof(double) * (x->fft_size));
            x->outs = (double *) fftw_malloc(sizeof(double) * (x->fft_size)*2);
            memset(x->in, '\0', x->fft_size * sizeof(double));
            memset(x->outs, '\0', x->fft_size * 2 * sizeof(double));
            if(x->mag_spec !=NULL) free(x->mag_spec);
            if(x->phase_spec !=NULL) free(x->phase_spec);
            x->mag_spec = malloc(sizeof(double)*(x->fft_size));     //these could be half as long
            x->phase_spec = malloc(sizeof(double)*(x->fft_size));
           
            if(x->ap0 !=NULL) free(x->ap0);
            x->ap0 = malloc(sizeof(double)*(x->fft_size / 2));
            
            if(x->ap1 !=NULL) free(x->ap1);
            x->ap1 = malloc(sizeof(double)*(x->fft_size / 2));
            
            if(x->ap2 !=NULL) free(x->ap2);
            x->ap2 = malloc(sizeof(double)*(x->fft_size / 2));
            
            if(x->ap3 !=NULL) free(x->ap3);
            x->ap3 = malloc(sizeof(double)*(x->fft_size / 2));
            
            if(x->ap4 !=NULL) free(x->ap4);
            x->ap4 = malloc(sizeof(double)*(x->fft_size / 2));
            
            //we also need to change the size of our cooked array:
            if(x->cooked !=NULL) free(x->cooked);
            x->cooked = malloc(sizeof(double) * (x->fft_size));
            
            //make a new plan. If one already exists, kill it and build a new one for the new fft size.
            //note, using the FFTW_MEASURE flag adds a beat, but optimizes for fast execution.
            //one could also use the FFTW_PATIENT flag here to really optimize at the expense an even longer wait.
            if(x->p !=NULL) fftw_destroy_plan(x->p);
            x->p = fftw_plan_dft_r2c_1d(x->fft_size, x->in, (fftw_complex *)x->outs, FFTW_MEASURE);
            
            

            post("jw_vindex: FFT size set to %d", n);
            
        } else {
            post("jw_vindex: FFT size must be a power of 2");
        }
    } else {
        post("jw_vindex: FFT size must be a positive integer");
    }
}

void jw_vindex_getvsize(t_jw_vindex *x)
{
    outlet_int(x->out, x->sample_vector_size);
    post("jw_vindex: Sample Vector Size is %d", x->sample_vector_size);
}

void jw_vindex_in1(t_jw_vindex *x, long n)
{
    if (n)
        x->l_chan = MAX(n, 1) - 1;
    else
        x->l_chan = 0;
}


void jw_vindex_dsp64(t_jw_vindex *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
    dsp_add64(dsp64, (t_object *)x, (t_perfroutine64)jw_vindex_perform64, 0, NULL);
}


// this lets us double-click on index~ to open up the buffer~ it references
void jw_vindex_dblclick(t_jw_vindex *x)
{
    buffer_view(buffer_ref_getobject(x->l_buffer_reference));
}

void jw_vindex_assist(t_jw_vindex *x, void *b, long m, long a, char *s)
{
    if (m == ASSIST_OUTLET)
        sprintf(s,"(signal) Sample Value at Index");
    else {
        switch (a) {
        case 0:    sprintf(s,"(signal) Sample Index");    break;
        case 1:    sprintf(s,"Audio Channel In buffer~");    break;
        case 3:    sprintf(s,"Dump Outlet");    break;
        }
    }
    
}

void jw_vindex_set_thresh(t_jw_vindex *x, double n)
{
    if(n<=0){
        x->thresh = n;
        post("jw_vindex: Threshold set to %f dB relative to peak value.", n);
    } else {
        post("jw_vindex: Threshold should be a negative number (dB relative to peak value)");
    }
}

//eventually, I'd like the bang method to find the first peak in the buffer and return the resonance at that point. For now, it outputs the last frame
void jw_vindex_bang(t_jw_vindex *x)
{
//    t_atom myList[3];
//    double theNumbers[3];
//    short i;
//
//    theNumbers[0] = 23.01;
//    theNumbers[1] = 12.02;
//    theNumbers[2] = 5.03;
//    for (i=0; i < 3; i++) {
//        atom_setfloat(myList+i,theNumbers[i]);
//    }
//    outlet_list(x->d_out, 0L, 3, &myList);
    jw_vindex_list_out(x, x->cooked, x->num_peaks*2);     //list the cooked frequencies out the
}

void jw_vindex_list_out(t_jw_vindex *x, double* a, long l)
{
    t_atom list[l];
    for(int i = 0; i<l; i++){
        atom_setfloat(list+i, a[i]);
    }
    outlet_list(x->d_out, 0L, l, list);
    
}

void *jw_vindex_new(t_symbol *s, long chan)
{
    t_jw_vindex *x = object_alloc(jw_vindex_class);
    dsp_setup((t_pxobject *)x, 1);
    intin((t_object *)x,1);
    x->out = outlet_new((t_object *)x, "int");  //right outlet
    x->f_out = outlet_new((t_object *)x, "float");
    x->d_out = outlet_new((t_object *)x, NULL);
    outlet_new((t_object *)x, "signal");        //left outlet
    jw_vindex_set(x, s);
    jw_vindex_in1(x,chan);
    x->sr = sys_getsr();                        //initially, adopt the system sample rate. We will reset later.
    post("jw_vindex: SR = %f", x->sr);
    x->fft_size = 1024;
    x->in = (double *) fftw_malloc(sizeof(double) * (x->fft_size));
    x->outs = (double *) fftw_malloc(sizeof(double) * (x->fft_size)*2);
    memset(x->in, '\0', x->fft_size * sizeof(double));  //initialize to zero
    memset(x->outs, '\0', x->fft_size * 2 * sizeof(double));
    x->p = fftw_plan_dft_r2c_1d(x->fft_size, x->in, (fftw_complex *)x->outs, FFTW_MEASURE);
    
    x->mag_spec = malloc(sizeof(double)*(x->fft_size));     //these could be half as long
    x->phase_spec = malloc(sizeof(double)*(x->fft_size));
    x->peaks = malloc(sizeof(long) * (x->fft_size / 2));
    x->cooked = malloc(sizeof(double) * (x->fft_size));
    x->analysis_points = malloc(sizeof(long) * 5);      //we are going to analyze just 5 points to extract decay rates
    x->ap0 = malloc(sizeof(double) * x->fft_size / 2);
    x->ap1 = malloc(sizeof(double) * x->fft_size / 2);
    x->ap2 = malloc(sizeof(double) * x->fft_size / 2);
    x->ap3 = malloc(sizeof(double) * x->fft_size / 2);
    x->ap4 = malloc(sizeof(double) * x->fft_size / 2);
    
    x->thresh = -32;
    x->num_peaks = 0;
    
    //test exponential fitting
    long *xVals = malloc(sizeof(long)*5);
    double *yVals = malloc(sizeof(double)*5);
    
    for(long i=1;i<6;i++){
        xVals[i-1]=i;
        yVals[i-1]= 3 * exp(-0.002 * (double)i)+0.005*random();
    }
    
    double *eftest = exp_fit(xVals,yVals, 5);
    post("jw_vindex: exp_fit test: y=A*exp(B*x) for vals (%d, %lf), (%d, %lf), (%d, %lf), (%d, %f), (%d, %f)",
         xVals[0], yVals[0],
         xVals[1], yVals[1],
         xVals[2], yVals[2],
         xVals[3], yVals[3],
         xVals[4], yVals[4]
         );
    post("A: %f", eftest[0]);
    post("B: %f", eftest[1]);   //works
    
    return (x);
}


void jw_vindex_free(t_jw_vindex *x)
{
    dsp_free((t_pxobject *)x);
    fftw_destroy_plan(x->p);
    if(x->in != NULL) fftw_free((char *)x->in);
    if(x->outs !=NULL) fftw_free((char *)x->outs);
    if(x->mag_spec !=NULL) free(x->mag_spec);
    if(x->phase_spec !=NULL) free(x->phase_spec);
    if(x->ap0 !=NULL) free(x->ap0);
    if(x->ap1 !=NULL) free(x->ap1);
    if(x->ap2 !=NULL) free(x->ap2);
    if(x->ap3 !=NULL) free(x->ap3);
    if(x->ap4 !=NULL) free(x->ap4);
    object_free(x->l_buffer_reference);
    
}


t_max_err jw_vindex_notify(t_jw_vindex *x, t_symbol *s, t_symbol *msg, void *sender, void *data)
{
    return buffer_ref_notify(x->l_buffer_reference, s, msg, sender, data);
}

//hann_window function
void hann_window(t_jw_vindex *x){
    int N = x->fft_size;
    for(int i=0;i<N;i++){
        double scale = pow(sin(PI*i/N),2);
        x->in[i] *= scale;
    }
}

//method to perform exponential fitting via least squares
double *exp_fit(long *xVals, double *yVals, long n)
{
    double *out = malloc(sizeof(double)*2);
    out[0] =0;
    out[1] = 0;
    double sum_Y=0;
    double sum_XY=0;
    double sum_X2Y=0;
    double sum_YlnY=0;
    double sum_XYlnY=0;
    
    for(int i=0;i<n;i++){
        double XY = (double)xVals[i] * yVals[i];
        double X2Y = ((double)xVals[i]) * XY;
        double YlnY = yVals[i] * log(yVals[i]);
        double XYlnY = (double)xVals[i] * YlnY;
        
        sum_Y += yVals[i];
        sum_XY += XY;
        sum_X2Y += X2Y;
        sum_YlnY += YlnY;
        sum_XYlnY += XYlnY;
    }
    post("jw_vindex:XY: %f, %f, %f, %f, %f", sum_Y, sum_XY, sum_X2Y, sum_YlnY, sum_XYlnY);
    
    
    double den = sum_Y * sum_X2Y - sum_XY * sum_XY;
    
    out[0] = exp((sum_X2Y * sum_YlnY - sum_XY * sum_XYlnY)/(den));
    out[1] = (sum_Y * sum_XYlnY - sum_XY * sum_YlnY)/(den);
    
    post("%f", out[0]);
    post("%f", out[1]);
    
    
    
    return out;
}
