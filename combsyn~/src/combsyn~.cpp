// combsyn~.cpp
//
// (c) 2012 www.carminecella.com
//

#include "Comb.h"
#include "maxmix.h"
#include "maxcpp5.h"
#include <string>
#include <vector>

const int MAX_VSIZE = 16384;

static t_symbol *sym_getname =  gensym((char*) "getname");

using namespace std;
using namespace combsyn;

float frand (float min, float max) {
	float f = 0;
	short r = (short) (rand ());
	f = fabs ((float) r / 32768);
	f *= (max - min);
	f += min;
	return f;
}

class CombSynMax : public MspCpp5<CombSynMax> {
public:
    //parameters    
    t_sample m_soundness;
    t_sample m_damping;
    t_sample m_norm;
    
    int m_N;
    std::vector<Comb<float>* > m_combs;
//    std::vector<Delay<float>* > m_delays;
    std::vector<float> m_freqs;
//    std::vector<float> m_pannings;
    float m_cnorm;
    
    float* obuff;
//    float* obuffS;
    
    CombSynMax(t_symbol * sym, long ac, t_atom * av) {
        setupIO (&CombSynMax::perform, 1, 1);

        //parameters initialisation        
        m_soundness = .999;
        m_damping = .1;
        m_norm = 1.;
        
        m_N = 1;

        if (ac == 1 && isNumber (av)) {
            m_N = getNumberInt (av);
        }
        else {
			object_error ((t_object*) this, (char*) "combsyn~::warning: missing number of voices; defaulting to 1");
            m_N = 1;
		}
        
        float freq = 110; // A3 by default
        float coeff = 1.67;
        float stot = 0;
        for (int i = 0; i < m_N; ++i) {
            int n = i + 1;
            float fcurr = freq * pow ((float) n, coeff);
            if (fcurr > sys_getsr() / 2) fcurr -= sys_getsr() / 2;
            float s = m_soundness; //frand (m_soundness - .05, m_soundness + .05); 
            stot += s;
            Comb<float>* cb = new Comb<float> (sys_getsr(), 1. / fcurr, s);
            m_combs.push_back (cb);
            m_freqs.push_back (fcurr);
            float da = m_damping; //frand (m_damping - .05, m_damping + .05); 
            cb->damping (da);
            
//            int d = (int) frand (11, 73);
//            m_delays.push_back(new Delay<float> (sys_getsr(), d, 0)); // stereo decorrelation	
//            float p = frand (.1, .9); // random panning
//            m_pannings.push_back (p);
//            
            post ("freq: %g, sonance: %g", fcurr, s);
        }
        
        obuff = new float[MAX_VSIZE];
//        obuffS = new float[MAX_VSIZE];

        m_cnorm = 1. / m_freqs.size ();
        //m_cnorm = 1. - (1. / m_freqs.size ());

		//if (m_cnorm < .1) m_cnorm = .1;
		
    }

    ~CombSynMax() {
        for (unsigned int i = 0; i < m_freqs.size (); ++i)  {
            delete m_combs[i];
        }
        delete [] obuff;
        //delete [] obuffS;
    }

    // optional method: gets called when the dsp chain is modified
    void dsp() {
    }
	void bang (long inlet) { 
        for (unsigned int j = 0; j < m_freqs.size (); ++j) {
            m_combs[j]->reset ();
        }
        post ("reset");
	}
    
   void tune (long inlet, t_symbol* s, long argc, t_atom* argv) { 
        if (argc > m_freqs.size () || argc < 1) {
            post ("voices submitted = %d (%s)", argc, s->s_name);
            for (int i = 0; i < argc; ++i)  {
                post ("argv %d = %g", i + 1, getNumberFloat (argv + i));
            }
           object_error ((t_object*) this, "combsyn~::error: invalid number of voices ");
           return;           
       }

       for (int i = 0; i < argc; ++i)  {
           delete m_combs[i];
           float freq = getNumberFloat (argv + i);
           m_combs[i] = new Comb<float> (sys_getsr(), 1. / freq, m_soundness);
           m_combs[i]->damping (m_damping);
           m_freqs[i] = freq;
           post ("new freq: %g", freq);
       }
   }

    void design (long inlet, t_symbol* s, long argc, t_atom* argv) { 
        if (argc != 2) {
            object_error ((t_object*) this, "combsyn~::error: synxat is '<float:f0>, <float:coeff>'");
            return;           
        }

        float freq = getNumberFloat (argv);
        float coeff = getNumberFloat (argv + 1);
        for (int i = 0; i < m_freqs.size (); ++i)  {
            delete m_combs[i];
            int n = i + 1;
            float fcurr = freq * pow ((float) n, coeff);
            if (fcurr > sys_getsr() / 2) fcurr -= sys_getsr() / 2;
            m_combs[i] = new Comb<float> (sys_getsr(), 1. / fcurr, m_soundness);
            m_combs[i]->damping (m_damping);
            m_freqs[i] = fcurr;
            post ("new freq: %g", fcurr);
        }
    }
    
    // signal processing CombSyn
    void perform (int vs, t_sample ** inputs, t_sample ** outputs) {
        t_sample *in = inputs[0];
        t_sample *out = outputs[0];
        
        memset (out, 0, sizeof (float) * vs);
        for (unsigned int j = 0; j < m_freqs.size (); ++j) {
            memset (obuff, 0, sizeof (float) * vs);
            m_combs[j]->process ((const float*)in, obuff, vs);
            
           // m_delays[j]->process (in, obuffS, sys_getblksize ());
            
//            for (int k = 0; k < vs; ++k) {
//                out[k * 2] = in[k]; //* m_pannings[j] * m_norm * m_cnorm;
//                out[k * 2 + 1] = in[k]; // * (1. - m_pannings[j]) * m_norm * m_cnorm;
//            }
            for (int k = 0; k < vs; ++k) {
                out[k] += obuff[k] * m_norm * m_cnorm;
            }

        }


    }
};

// attributes
t_max_err getAttr (CombSynMax *self, t_object *attr, long* ac, t_atom** av) {
    if ((*ac) == 0 || (*av) == NULL) {
        //otherwise allocate memory
        *ac = 1;

        if (!(*av = (t_atom *)getbytes (sizeof (t_atom) * (*ac)))) {
            *ac = 0;
            return MAX_ERR_OUT_OF_MEM;
        }
    }

    string attrname = ((t_symbol *)object_method ((t_object *)attr, sym_getname))->s_name;

    if (attrname.compare ("soundness") == 0) {
        atom_setfloat(*av, self->m_soundness);
        return MAX_ERR_NONE;
    }
    
    if (attrname.compare ("damping") == 0) {
        atom_setfloat(*av, self->m_damping);
        return MAX_ERR_NONE;
    }    
    
    if (attrname.compare ("normalization") == 0) {
        atom_setfloat(*av, self->m_norm);
        return MAX_ERR_NONE;
    }

    return MAX_ERR_NONE;
}


t_max_err setAttr (CombSynMax *self, void *attr, long ac, t_atom *av) {
    if (!(ac>0 && isNumber(av))) return MAX_ERR_NONE;

    t_symbol *attrsym = (t_symbol *)object_method((t_object *)attr, sym_getname);
    string attrname = attrsym->s_name;

    if (attrname.compare("soundness") == 0) {
        self->m_soundness = getNumberFloat(av);
    }
    
    if (attrname.compare("damping") == 0) {
        self->m_damping = getNumberFloat(av);
    }
    
    if (attrname.compare("normalization") == 0) {
        self->m_norm = getNumberFloat(av);
    }

    float stot = 0;
    for (unsigned int i = 0; i < self->m_freqs.size (); ++i) {
//        float s = frand (self->m_soundness - .05, self->m_soundness + 0.05);
//        float da = frand (self->m_damping - .05, self->m_damping + 0.05);
//        if (s > .9999) s = .9999;
//        if (da > .9999) da = .9999;
        stot += self->m_soundness;
        self->m_combs[i]->feedback (self->m_soundness);
        self->m_combs[i]->damping (self->m_damping);
    }
//    
//    post ("soundness: %g, damping: %g", self->m_soundness, self->m_damping);
    
//    self->m_cnorm = 1. - (stot / self->m_freqs.size ());
//    if (self->m_cnorm < .1) self->m_cnorm = .1;
    
    return MAX_ERR_NONE;
}

extern "C" int main(void) {
    // create a class with the given name:
    CombSynMax::makeMaxClass("combsyn~");
    t_class* c = (t_class *)CombSynMax::m_class;

    CLASS_ATTR_FLOAT(c, (char*) "soundness", 0, CombSynMax, m_soundness);
    CLASS_ATTR_ACCESSORS(c, (char*) "soundness", (method)getAttr, (method)setAttr);
    CLASS_ATTR_FILTER_MIN(c, (char*) "soundness", 0);
    CLASS_ATTR_SAVE(c, (char*) "soundness", .999);
    
    CLASS_ATTR_FLOAT(c, (char*) "damping", 0, CombSynMax, m_damping);
    CLASS_ATTR_ACCESSORS(c, (char*) "damping", (method)getAttr, (method)setAttr);
    CLASS_ATTR_FILTER_MIN(c, (char*) "damping", 0.001);
    CLASS_ATTR_SAVE(c, (char*) "damping", .1);
    
    CLASS_ATTR_FLOAT(c, (char*) "normalization", 0, CombSynMax, m_norm);
    CLASS_ATTR_ACCESSORS(c, (char*) "normalization", (method)getAttr, (method)setAttr);
    CLASS_ATTR_FILTER_MIN(c, (char*) "normalization", 0);
    CLASS_ATTR_SAVE(c, (char*) "normalization", .1);
        
    REGISTER_METHOD(CombSynMax, bang);
	REGISTER_METHOD_GIMME(CombSynMax, tune);
   	REGISTER_METHOD_GIMME(CombSynMax, design);

}


// eof

