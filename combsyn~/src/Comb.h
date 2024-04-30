// Comb.h
// 


#ifndef COMB_H
#define COMB_H 

#include "Delay.h"
namespace combsyn  {
template <typename T>
//! Comb filter: for each buffer returns a buffer of same size
class Comb : public Delay <T> {
public:
    typedef Delay <T> Base;
    Comb (double sr, double dur, double feedback) : 
        Delay <T> (sr, dur, feedback), 
        m_damp (0), m_udamp (1), m_lp (0) {}
    Comb (double sr, int samples, double feedback) : 
        Delay <T> (sr, samples, feedback), 
        m_damp (0), m_udamp (1), m_lp (0) {}		
    void damping (double d, bool reset = false) {
        m_damp = d;
        m_udamp = 1. - m_damp;
        if (reset) {
            m_lp = 0;
        }
    }
    virtual T* process (const T* input, T* output, long size) {
        for (long i = 0; i < size; ++i) {
            output[i] = Base::m_delay[Base::m_ptr];
            //undenormalise (output[i]);
            m_lp = (output[i] * m_udamp) + (m_lp * m_damp);
            //undenormalise (m_lp);
            Base::m_delay[Base::m_ptr] = input[i] + (Base::m_feedback * m_lp);
            Base::m_ptr++;
            if (Base::m_ptr < 0) Base::m_ptr = 0;
            Base::m_ptr %= Base::m_samples;	
        }
        
        return output;
    }
private:
    double m_lp;
    double m_damp;
    double m_udamp;
};
}
#endif	// COMB_H 

// EOF
