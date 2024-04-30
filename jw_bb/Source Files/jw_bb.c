//
//  jw_bb.c
//
//
//  Created by Jeremy Wagner on 4/24/24.
//

#include "jw_bb.h"
/**
    @file
    jw_bb - a bollinger band generator
    Jeremy Wagner

    @ingroup    examples
*/
#include "ext.h"                            // standard Max include, always required
#include "ext_obex.h"                        // required for new style Max object

#define RING_BUFFER_SIZE 256
////////////////////////// object struct
typedef struct _jw_bb
{
    t_object                    ob;             // the object itself (must be first)
    void*                       out_ubb;        //left outlet; upper bollinger band
    void*                       out_lbb;        //middle outlet; lower bollinger band
    void*                       out_running_avg;    //right outlet; running average
    t_atom_long                 num_samples;    //number of samples in the running average
    double                      num_sd;         //number of standard deviations above/below
    double*                     ring_buffer;    //buffer of values
    t_atom_long                 cursor;         //buffer write position
    double                      mean;           //running avg
    double                      stdDev;       //stdDev
    double                      ubb;          //upper Bollinger Band
    double                      lbb;          //lower Bollinger Band
} t_jw_bb;

///////////////////////// function prototypes
//// standard set
void *jw_bb_new(t_symbol *s, long argc, t_atom *argv);
void jw_bb_free(t_jw_bb *x);
void jw_bb_assist(t_jw_bb *x, void *b, long m, long a, char *s);
void jw_bb_int(t_jw_bb *x, long n);
void jw_bb_float(t_jw_bb *x, double f);
//void jw_bb_list(t_jw_bb *x, t_symbol *msg, long argc, t_atom *argv);
void jw_bb_set_num(t_jw_bb *x, long n);
void jw_bb_num_sd(t_jw_bb *x, double f);
void jw_bb_clear(t_jw_bb *x);


//////////////////////// global class pointer variable
void *jw_bb_class;


void ext_main(void *r)
{
    t_class *c;

    c = class_new("jw_bb", (method)jw_bb_new, (method)jw_bb_free, (long)sizeof(t_jw_bb),
                  0L /* leave NULL!! */, A_GIMME, 0);
    
    /* you CAN'T call this from the patcher */
    class_addmethod(c, (method)jw_bb_assist,            "assist",        A_CANT, 0);
    class_addmethod(c, (method)jw_bb_int,               "int",           A_LONG, 0);
    class_addmethod(c, (method)jw_bb_float,             "float",         A_FLOAT, 0);
//    class_addmethod(c, (method)jw_bb_list,              "list",          A_CANT, 0);
    class_addmethod(c, (method)jw_bb_set_num,           "num_samples",   A_LONG, 0);
    class_addmethod(c, (method)jw_bb_num_sd,            "num_sd",        A_FLOAT, 0);
    class_addmethod(c, (method)jw_bb_clear,             "clear",         A_DEFSYM, 0);

    class_register(CLASS_BOX, c); /* CLASS_NOBOX */
    jw_bb_class = c;

    post("I am the jw_bb object");
}

static void NewInput(t_jw_bb *x, double f)
{
    x->ring_buffer[x->cursor] = f;
    ++(x->cursor);
    x->cursor %= x->num_samples;
}

void bollingerBands(t_jw_bb *x)
{
    double mu = 0;
    for(int i=0;i<x->num_samples; i++){mu += x->ring_buffer[i];}
    x->mean = mu / x->num_samples;
    double sum = 0;
    for(int i=0;i<x->num_samples; i++){
        sum += (x->ring_buffer[i] - x->mean) * (x->ring_buffer[i] - x->mean);
    }
    x->stdDev = sqrt(sum / x->num_samples);
    x->ubb = x->mean + x->num_sd * x->stdDev;
    x->lbb = x->mean - x->num_sd * x->stdDev;
}

void jw_bb_set_num(t_jw_bb *x, long n)
{
    if(n>RING_BUFFER_SIZE){
        post("jw_bb: num_samples limited to %d. Setting to %d", RING_BUFFER_SIZE, RING_BUFFER_SIZE);
        x->num_samples = RING_BUFFER_SIZE;
    }else{
        x->num_samples = n;
    }
    jw_bb_clear(x);
    x->cursor = 0;
}

void jw_bb_num_sd(t_jw_bb *x, double f)
{
    x->num_sd = fabs(f);
}

void jw_bb_clear(t_jw_bb *x)
{
    for(int i=0; i<RING_BUFFER_SIZE;i++) x->ring_buffer[i] = 0;
}

void jw_bb_float(t_jw_bb *x, double f)
{
    NewInput(x, f);
    bollingerBands(x);
    outlet_float(x->out_lbb, x->lbb);
    outlet_float(x->out_ubb, x->ubb);
    outlet_float(x->out_running_avg, x->mean);
}

void jw_bb_int(t_jw_bb *x, long n)
{
    NewInput(x, (double)n);
    bollingerBands(x);
    outlet_float(x->out_lbb, x->lbb);
    outlet_float(x->out_ubb, x->ubb);
    outlet_float(x->out_running_avg, x->mean);


}

void jw_bb_assist(t_jw_bb *x, void *b, long m, long a, char *s)
{
    if (m == ASSIST_INLET) { // inlet
        sprintf(s, "I am inlet %ld", a);
    }
    else {    // outlet
        sprintf(s, "I am outlet %ld", a);
    }
}

void jw_bb_free(t_jw_bb *x)
{
    if(x->ring_buffer !=NULL) free(x->ring_buffer);
}


void *jw_bb_new(t_symbol *s, long argc, t_atom *argv)
{
    t_jw_bb *x = NULL;
    long i;

    if ((x = (t_jw_bb *)object_alloc(jw_bb_class))) {
        object_post((t_object *)x, "a new %s object was instantiated: %p", s->s_name, x);
        object_post((t_object *)x, "it has %ld arguments", argc);

        for (i = 0; i < argc; i++) {
            if ((argv + i)->a_type == A_LONG) {
                object_post((t_object *)x, "arg %ld: long (%ld)", i, atom_getlong(argv+i));
                x->num_samples=atom_getlong(argv+i);
            } else if ((argv + i)->a_type == A_FLOAT) {
                object_post((t_object *)x, "arg %ld: float (%f)", i, atom_getfloat(argv+i));
                x->num_sd=atom_getlong(argv+i);
            } else {
                object_error((t_object *)x, "forbidden argument");
                x->num_samples = RING_BUFFER_SIZE;
            }
        }
    }
    
    x->ring_buffer = (double*)calloc(x->num_samples, sizeof(double));
    x->mean = 0.0;
    x->lbb=0.0;
    x->ubb=0.0;
    
    
    x->out_running_avg = outlet_new((t_object *)x, "float");  //right outlet; running average
    x->out_lbb = outlet_new((t_object *)x, "float");          //middle outlet; lower bollinger band
    x->out_ubb = outlet_new((t_object *)x, "float");          //left outlet; upper bollinger band
    
    
    return (x);
}
