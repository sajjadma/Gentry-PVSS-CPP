/* my_encryption_type2.cpp
 *
 * Copyright (C) 2025
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#include <my_implementation.hpp>

namespace MyEncryption {
    // TODO: Set Parameters Properly
    void EncryptionType2::setup(const MyFramework::Encryption::Params *params, const long securityParameter,
                                const NTL::ZZ plainBound) {
        const auto myParams = (MyParams *) params;
        myParams->k = NTL::NumBits(plainBound);
        myParams->k += myParams->k % 2;
        myParams->module = NTL::power2_ZZ(myParams->k);
        myParams->l = myParams->k; // TODO: یا هر مقدار مناسب دیگر
        myParams->m = myParams->l * myParams->k;
        myParams->d = 2; // TODO: یا هر مقدار مناسب دیگر
        myParams->randomBound = NTL::power2_ZZ(myParams->k / 2 - 2);
        myParams->plainBound = myParams->module;
        myParams->coefficientBound = myParams->module;
        myParams->plainSize = 1;
        myParams->randomSize = myParams->m * myParams->d + myParams->l + myParams->m + myParams->d;
        myParams->cipherSize = myParams->l * myParams->d + myParams->m + myParams->d;
    }
}
