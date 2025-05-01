#ifndef PVSSSCHEMETYPE1_HPP
#define PVSSSCHEMETYPE1_HPP
/* pvssSchemeType1.hpp
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

#include <vector>
#include "NTL/mat_ZZ_pE.h"

namespace NewPVSSScheme::PVSSType1 {
    struct Params {
        long numberOfParties;
        long threshold;
        NTL::ZZ p;
        NTL::ZZ q;
        long k;
        NTL::ZZ_pX f;
        NTL::ZZ bound;
        NTL::ZZ_pE a;
        NTL::vec_ZZ_pE v;
        NTL::vec_ZZ_pE h;
        NTL::mat_ZZ_pE A;
        NTL::vec_ZZ_pE t;
        NTL::Mat<NTL::vec_ZZ_pE> u;
    };

    struct PrivateKey {
        NTL::ZZ_pE a;
        NTL::ZZ_pE s;
    };

    struct PublicKey {
        NTL::ZZ_pE a;
        NTL::ZZ_pE b;
    };

    struct KeyProof {
    };

    struct KeyPair {
        PublicKey publicKey;
        PrivateKey privateKey;
        KeyProof proof;
    };

    struct Cipher {
        NTL::ZZ_pE u;
        NTL::ZZ_pE v;
    };

    struct Commitment {
        NTL::ZZ_pE c;
    };

    struct OpeningProof {
        NTL::vec_ZZ output;
    };

    struct DistributionProof {
        std::vector<Cipher> encryptedShares;
        Commitment commitment;
        OpeningProof proof;
    };

    struct DecryptionProof {
        NTL::ZZ decryptedShare;
    };

    void _generateTrapdoor(NTL::mat_ZZ_pE &A, NTL::mat_ZZ_pE &trapdoor, long n, long m,
                           const NTL::ZZ &q, const NTL::ZZ_pX &f, const NTL::ZZ &bound);

    void _preSample(NTL::vec_ZZ_pE &x, const NTL::mat_ZZ_pE &trapdoor, const NTL::mat_ZZ_pE &A, const NTL::vec_ZZ_pE &b,
                    const NTL::ZZ &q, const NTL::ZZ_pX &f, const NTL::ZZ &bound);

    Params setup(long securityParameter, long numberOfParties, long threshold);

    KeyPair generateKey(const Params &params, long index);

    bool verifyKey(const Params &params, const PublicKey &publicKey, const KeyProof &proof);

    DistributionProof distribute(const Params &params, const std::vector<PublicKey> &publicKeys,
                                 const NTL::ZZ &secret);

    bool verifyDistribution(const Params &params, const std::vector<PublicKey> &publicKeys,
                            const DistributionProof &proof);

    DecryptionProof decryptShare(const Params &params, const PublicKey &publicKey,
                                 const PrivateKey &privateKey, const Cipher &encryptedShare);

    bool verifyDecryption(const Params &params, const PublicKey &publicKey,
                          const Cipher &encryptedShare, const DecryptionProof &proof);

    NTL::ZZ reconstruct(const Params &params, const std::vector<NTL::ZZ> &decryptedShares);
}

#endif //PVSSSCHEMETYPE1_HPP
