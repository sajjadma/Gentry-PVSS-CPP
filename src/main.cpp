/* main.cpp - a "main" file, just a debugging tool
 * 
 * Copyright (C) 2021, LWE-PVSS
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
 
// This file is just a convenience, a handy tool that lets us run
// small porgrams without having to use the awkward ctest syntax.
#include <cassert>
#include <cmath>
#include <random>
#include <chrono>
#include <my_implementation.hpp>
#include <pvss_framework.hpp>
#include <string>
#include <sys/time.h>
#include <sys/resource.h>
using namespace std;

#include <NTL/version.h>
#include "regevEnc.hpp"
#include "regevProofs.hpp"

using namespace ALGEBRA;
using namespace REGEVENC;

int main(int argc, char** argv) {
    // std::cout << "- Found GMP version "<<__GNU_MP__ <<std::endl;
    // MyVectorCommitment::VectorCommitmentType1 vc;
    // MyFramework::VC::Params *vcParams = new MyVectorCommitment::VectorCommitmentType1::MyParams();
    // NTL::ZZ prime = NTL::conv<NTL::ZZ>("27742317777372353535851937790883648493");
    // NTL::ZZ bound = NTL::power2_ZZ(8);
    // vc.setup(vcParams, 1, 2, 4, 8, prime, bound, prime);


    auto start2 = chrono::steady_clock::now();
    cout << "PVSS starts" << endl;

    int n = 5;
    MyFramework::PVSS pvss(
        make_unique<MyEncryption::EncryptionType1>(),
        make_unique<MyVectorCommitment::VectorCommitmentType2>()
    );
    MyFramework::Params params;
    params.encryptionParams = new MyEncryption::EncryptionType1::MyParams();
    params.vcParams = new MyVectorCommitment::VectorCommitmentType2::MyParams();
    pvss.setup(params, 1, n, n/2 + 1);


    auto start1 = chrono::steady_clock::now();
    cout << "PVSS KeyGen All starts" << endl;

    MyFramework::Encryption::KeyPair keyPair[n];
    vector<MyFramework::Encryption::PublicKey *> pks;
    for (auto &key: keyPair) {
        key.privateKey = new MyEncryption::EncryptionType1::MyPrivateKey();
        key.publicKey = new MyEncryption::EncryptionType1::MyPublicKey();
        key.proof = new MyEncryption::EncryptionType1::MyKeyProof();
        pvss.generateKey(key, params);
        pks.push_back(key.publicKey);
    }

    auto end1 = chrono::steady_clock::now();
    auto ticks1 = chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();
    cout << "PVSS KeyGen All ends. time: " << ticks1 << "ms" << endl;

    cout << "verifyKey: " << pvss.verifyKey(params, keyPair[n-1].publicKey, keyPair[n-1].proof) << endl;

    MyFramework::DistributionProof proof;
    proof.commitment = new MyVectorCommitment::VectorCommitmentType2::MyCommitment();
    proof.proof = new MyVectorCommitment::VectorCommitmentType2::MyOpeningProof();
    pvss.distribute(proof, params, pks, NTL::to_ZZ(5));

    cout << "verifyDistribution: " << pvss.verifyDistribution(params, pks, proof) << endl;


    start1 = chrono::steady_clock::now();
    cout << "PVSS Dec All starts" << endl;

    MyFramework::DecryptionProof decrypt;
    decrypt.proof = new MyEncryption::EncryptionType1::MyDecryptionProof();
    vector<NTL::ZZ> sh;
    for (int i = 0; i < n; i++) {
        pvss.decryptShare(decrypt,params,keyPair[i].publicKey,keyPair[i].privateKey,proof.encryptedShares[i]);
        sh.push_back(decrypt.decryptedShare);
    }

    end1 = chrono::steady_clock::now();
    ticks1 = chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();
    cout << "PVSS Dec All ends. time: " << ticks1 << "ms" << endl;

    cout << "verifyDecryption: " << pvss.verifyDecryption(params, pks[n-1], proof.encryptedShares[n-1], decrypt) << endl;


    pvss.reconstruct(decrypt.decryptedShare, params, sh);

    cout << decrypt.decryptedShare;

    auto end2 = chrono::steady_clock::now();
    auto ticks2 = chrono::duration_cast<chrono::milliseconds>(end2 - start2).count();
    cout << "PVSS ends. time: " << ticks2 << "ms" << endl;

    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    std::cout << " max mem: " << ru.ru_maxrss << " kilobytes\n";


    return 0;
    // std::cout << "- Found NTL version "<<NTL_VERSION <<std::endl;
    // std::cout << "- Found Sodium version "<<SODIUM_VERSION_STRING<<std::endl;

    int nParties = 512;
    if (argc > 1) {
        nParties = std::stoi(argv[1]);
    }
    if (nParties < 32 || nParties > 4096)
        nParties = 512;
    std::cout << "nParties="<<nParties << std::endl;

    // The dimensions of the the CRX is k-by-m, but note that this is a
    // matrix over GF(p^2) so the lattice dimensions we get is twice that
    KeyParams kp(nParties);
    // kp.k=64;
    GlobalKey gpk("testContext", kp);
    std::cout <<"{ kay:"<<gpk.kay<<", enn:"<<gpk.enn
      <<", sigmaEnc1:"<<gpk.sigmaEnc1<<", sigmaEnc2:"<<gpk.sigmaEnc2<<" }\n";

    TernaryEMatrix::init();
    MerlinRegev mer;
    PedersenContext ped;
    SharingParams ssp(interval(1,gpk.enn+1), gpk.tee);
    VerifierData vd(gpk, ped, mer, ssp);
    ProverData pd(vd);

    // Generate/verify the proofs by the second party (idx=1)
    int partyIdx = 1;

    // Key generation for all the parties
    std::vector<ALGEBRA::EVector> kgNoise(gpk.enn);
    std::vector<ALGEBRA::EVector> sk(gpk.enn);
    std::vector<ALGEBRA::EVector> pk(gpk.enn);
    auto start = chrono::steady_clock::now();
    crsTicks = 0;
    for (int i=0; i<gpk.enn; i++) {
        std::tie(sk[i],pk[i]) = gpk.genKeys(&kgNoise[i]);
        gpk.addPK(pk[i]);
    }
    gpk.setKeyHash();
    auto end = chrono::steady_clock::now();
    auto ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    std::cout <<gpk.enn<<" keyGens in "<<ticks<<" milliseconds, avg="<<(ticks/double(gpk.enn))
        << " ("<< (crsTicks/double(gpk.enn)) << " for s x A)\n";

    // encryption
    std::vector<ALGEBRA::SVector> ptxt1(gpk.enn);
    std::vector<GlobalKey::CtxtPair> ctxt1(gpk.enn);
    // secret sharing of a random value , the secret itself is sshr[0]
    ALGEBRA::SVector sshr;
    ssp.randomSharing(sshr);
    for (int i=0; i<gpk.enn; i++) {
        resize(ptxt1[i], gpk.enn);
        for (int j=0; j<gpk.enn; j++) ptxt1[i][j] = sshr[i+1];
    }
    start = chrono::steady_clock::now();
    crsTicks = 0;
    for (int i=0; i<gpk.enn; i++) {
        ctxt1[i] = gpk.encrypt(ptxt1[i]);
    }
    end = chrono::steady_clock::now();
    ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    std::cout <<gpk.enn<<" encryptions in "<<ticks<<" milliseconds, avg="<<(ticks/double(gpk.enn)) 
        << " ("<< (crsTicks/double(gpk.enn)) << " for A x r)\n";

    // decryption at party #1
    ALGEBRA::SVector ptxt2;    resize(ptxt2, gpk.tee);
    ALGEBRA::EVector decNoise; resize(decNoise, gpk.tee);
    start = chrono::steady_clock::now();
    for (int i=0; i<gpk.tee; i++) { // decrypt 2nd entry in i'th ctxt
        ptxt2[i] = gpk.decrypt(sk[partyIdx], partyIdx, ctxt1[i], &(decNoise[i]));
    }
    end = chrono::steady_clock::now();
    ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    std::cout << gpk.tee << " decryptions in "<<ticks<< " milliseconds, avg="
        << (ticks/double(gpk.tee)) << std::endl;

    for (int i=0; i<gpk.tee; i++) { // decrypt 2nd entry in i'th ctxt
        if (ptxt2[i] != ptxt1[i][partyIdx])
            std::cout << "decryption error in "<<i<<"th ciphertext\n";
    }

    // re-encryption at party #1
    ALGEBRA::SVector ptxt3;
    resize(ptxt3, gpk.enn);
    for (int j=0; j<gpk.enn; j++) ptxt3[j] = sshr[j+1];
    ALGEBRA::EVector encRnd;
    REGEVENC::GlobalKey::CtxtPair eNoise;
    auto ctxt2 = gpk.encrypt(ptxt3, encRnd, eNoise);

    // Copy the first t ciphertexts into a k x t matrix and another t-vector
    EMatrix ctxtMat;
    resize(ctxtMat, gpk.kay, gpk.tee);
    EVector ctxtVec;
    resize(ctxtVec, gpk.tee);
    for (int i=0; i<gpk.tee; i++) {
        for (int j=0; j<gpk.kay; j++)
            ctxtMat[j][i] = ctxt1[i].first[j];
        ctxtVec[i] = ctxt1[i].second[partyIdx];
    }

    // prepare for proof, commit to the secret key
    DLPROOFS::Point::counter = 0;
    DLPROOFS::Point::timer = 0;
    int origSize = sk[partyIdx].length(); 
   
    vd.sk1Com = commit(sk[partyIdx], vd.sk1Idx, vd.Gs, pd.sk1Rnd);

    start = chrono::steady_clock::now();
    SVector lagrange = vd.sp->lagrangeCoeffs(interval(1,gpk.tee+1));

    crsTicks = 0;
    proveDecryption(pd, ctxtMat, ctxtVec, ptxt2, sk[partyIdx], decNoise);
    proveEncryption(pd, ctxt2.first, ctxt2.second, ptxt3, encRnd, eNoise.first, eNoise.second);
    proveKeyGen(pd, partyIdx, sk[partyIdx], kgNoise[partyIdx]);
    proveReShare(pd, lagrange, ptxt2, ptxt3);
    proveSmallness(pd);

    end = chrono::steady_clock::now();
    ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    std::cout << "preparing to prove and committing in "<<ticks<< " milliseconds, "
        << DLPROOFS::Point::counter << " exponentiations in "
        << ((500+DLPROOFS::Point::timer)/1000) << " milliseconds"
        << " ("<< crsTicks << " for xR x A)\n";

    // aggregate the constraints and flatten everything before proving
    DLPROOFS::Point::counter = 0;
    DLPROOFS::Point::timer = 0;
    std::cout<<"aggregting constraints\n";
    start = chrono::steady_clock::now();
    ReadyToProve rtp;
    rtp.aggregateProver(pd);

    // Make copies of the Merlin transcripts and specialize them
    // for the final constraints before proving/verifying them
    auto merLin = *vd.mer;
    merLin.processConstraint("linear", rtp.linCnstr);

    auto merQuad = *vd.mer;
    merQuad.processConstraint("quadratic", rtp.quadCnstr);

    // Flatten the statements, this relases the memory of the constraints
    // (hence the Merlin processing above must be done before doing this).
    std::cout<<"flatenning constraints\n";
    rtp.flattenLinPrv(pd);
    rtp.flattenQuadPrv(pd);

    ReadyToVerify rtv = rtp; // a copy without the secret variables

    // prove and verify the linear statement
    auto merLinVer = merLin; // another copy for verification
    DLPROOFS::LinPfTranscript pfL("Linear");
    pfL.C = rtp.linCom;

    end = chrono::steady_clock::now();
    ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    std::cout << "aggregating constaints in "<<ticks<< " milliseconds, "
        << DLPROOFS::Point::counter << " exponentiations in "
        << ((500+DLPROOFS::Point::timer)/1000) << " milliseconds\n";

    // The actual proof
    DLPROOFS::Point::counter = 0;
    DLPROOFS::Point::timer = 0;
    start = chrono::steady_clock::now();
    DLPROOFS::proveLinear(pfL, rtp.lComRnd, merLin, rtp.linWtns.data(),
            rtp.linStmnt.data(), rtp.linGs.data(), rtp.linGs.size());
    end = chrono::steady_clock::now();
    ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    std::cout << "proving linear in "<<ticks<< " milliseconds, "
        << DLPROOFS::Point::counter << " exponentiations in "
        << ((500+DLPROOFS::Point::timer)/1000) << " milliseconds\n";

    DLPROOFS::Point::counter = 0;
    DLPROOFS::Point::timer = 0;
    start = chrono::steady_clock::now();
    if (!DLPROOFS::verifyLinear(pfL, rtv.linStmnt.data(), rtv.linGs.data(),
                      rtv.linGs.size(), rtv.linCnstr.equalsTo, merLinVer))
        std::cout << "failed linear verification\n";
    end = chrono::steady_clock::now();
    ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    std::cout << "verifying linear in "<<ticks<< " milliseconds, "
        << DLPROOFS::Point::counter << " exponentiations in "
        << ((500+DLPROOFS::Point::timer)/1000) << " milliseconds\n";

    // prove and verify the quadratic statement
    DLPROOFS::Point::counter = 0;
    DLPROOFS::Point::timer = 0;
    auto merQuadVer = merQuad; // another copy for verification
    DLPROOFS::QuadPfTranscript pfQ("Quadratic");
    pfQ.C = rtp.quadCom;
     // The actual proof
    start = chrono::steady_clock::now();
    DLPROOFS::proveQuadratic(pfQ, rtp.qComRnd, merQuad, rtp.quadGs.data(),
                rtp.quadWtnsG.data(), rtp.quadHs.data(), rtp.quadWtnsH.data(),
                rtp.quadGs.size());
    end = chrono::steady_clock::now();
    ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    std::cout << "proving quadratic in "<<ticks<< " milliseconds, "
        << DLPROOFS::Point::counter << " exponentiations in "
        << ((500+DLPROOFS::Point::timer)/1000) << " milliseconds\n";

    // The actual verification
    DLPROOFS::Point::counter = 0;
    DLPROOFS::Point::timer = 0;
    start = chrono::steady_clock::now();
    if (!DLPROOFS::verifyQuadratic(pfQ, rtv.quadGs.data(), rtv.quadHs.data(),
                        rtp.quadGs.size(), rtv.quadCnstr.equalsTo, merQuadVer,
                        rtv.offstG.data(), rtv.offstH.data()))
        std::cout << "failed quadratic verification\n";
    end = chrono::steady_clock::now();
    ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    std::cout << "verifying quadratic in "<<ticks<< " milliseconds, "
        << DLPROOFS::Point::counter << " exponentiations in "
        << ((500+DLPROOFS::Point::timer)/1000) << " milliseconds\n";

    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    std::cout << " max mem: " << ru.ru_maxrss << " kilobytes\n";

    return 0;
}

#if 0
#include <vector>
#include <iostream>
#include "bulletproof.hpp"

int main(int, char**) {
    constexpr size_t pfSize = 13;

    // build a constraint: sum_i ai*bi = b = \sum_bi^2
    DLPROOFS::LinConstraint cnstrL;
    for (size_t i=0; i<pfSize; i++) {
        CRV25519::Scalar& x = cnstrL.terms[i+1].setInteger(i+1);
        cnstrL.equalsTo += x * x;
    }
    DLPROOFS::PtxtVec& xes = cnstrL.terms;
    DLPROOFS::LinPfTranscript pfL = proveLinear("blah", cnstrL, xes);
    std::cout << "linear: "<<verifyLinear(cnstrL, pfL) <<std::endl;

    DLPROOFS::QuadConstraint cnstrQ;
    for (auto& x : xes) {
        cnstrQ.indexes.insert(cnstrQ.indexes.end(), x.first);
        cnstrQ.equalsTo += x.second * x.second;
    }    
    DLPROOFS::PtxtVec& ys = cnstrL.terms;
    DLPROOFS::QuadPfTranscript pfQ = proveQuadratic("blah", cnstrQ, xes, ys);
    std::cout << "quadratic: "<<verifyQuadratic(cnstrQ, pfQ) <<std::endl;

    std::set<size_t> indexes;
    for (auto& elem: xes) // elem = {idx:scalar}
        indexes.insert(indexes.end(), elem.first);

    auto [normSq, prNS] = DLPROOFS::proveNormSquared("blah", xes);
    std::cout << "norm: "<<verifyNormSquared(indexes,normSq,prNS)<<std::endl;

    return 0;
}
#endif
