// Theres a few classes which I want to prohibit being created on the stack
// We use this struct to lock out the public contructors
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef KEY_HPP
#define KEY_HPP

#include <memory>

namespace analyticRT
{
    class raw_isobar;

    class key
    {
        private:

        // Private constructor only accessable via friend methods below
        key(){};

        template<class A>
        friend std::shared_ptr<raw_isobar> new_isobar();
        template<class A, class B>
        friend std::shared_ptr<raw_isobar> new_isobar(B);
        template<class A, class B, class C>
        friend std::shared_ptr<raw_isobar> new_isobar(B, C);
        template<class A, class B, class C, class D>
        friend std::shared_ptr<raw_isobar> new_isobar(B, C, D);
        template<class A, class B, class C, class D, class E>
        friend std::shared_ptr<raw_isobar> new_isobar(B, C, D, E);
    };
};

#endif