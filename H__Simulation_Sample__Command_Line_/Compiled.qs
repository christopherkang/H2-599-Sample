// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Samples.H2Simulation {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;

    // extension of the original H2Simulation file with compiled versions

    // Helper method to visualize the EvolutionGenerator
    operation ExploreGenerator(generator : EvolutionGenerator) : Unit {
        // let generator = H2EvolutionGenerator(2);
        let (evolSet, genSys) = generator!;
        let (maxNum, genIndexGen) = genSys!;

        for (index in 0..maxNum - 1) {
            Message($"{genIndexGen(index)}");
        }
    }

    /// # Summary
    /// We finish by using the Robust Phase Estimation algorithm
    /// of Kimmel, Low and Yoder.
    operation H2EstimateCompiledEnergyRPE(idxBondLength : Int, nBitsPrecision : Int, trotterStepSize : Double) : Double {
        return H2EstimateCompiledEnergy(idxBondLength, trotterStepSize, RobustPhaseEstimation(nBitsPrecision, _, _));
    }  

    // compiled energy
    // Note - the trotterStepsize has been forced to 1, as dictated by the protocol
    operation H2EstimateCompiledEnergy(idxBondLength : Int, trotterStepSize : Double, phaseEstAlgorithm : ((DiscreteOracle, Qubit[]) => Double)) : Double {

        let nQubits = 2;
        let trotterOrder = 1;

        //  constant - number of Hamiltonian terms to use
        let MAXSTEPS = 1;

        let newEvolutionGenerator = ConvertGenerator(MAXSTEPS, H2EvolutionGenerator(idxBondLength));
        // ExploreGenerator(newEvolutionGenerator);

        let trotterStep = H2CompiledTrotterStep(idxBondLength, trotterOrder, 1.0, _, newEvolutionGenerator);
        let estPhase = EstimateEnergy(nQubits, H2StatePrep, trotterStep, phaseEstAlgorithm);
        return estPhase / 1.0 + H2IdentityCoeff(idxBondLength);
    }

    /// # Summary
    /// We now provide Canon's Hamiltonian simulation
    /// functions with the above representation to automatically
    /// decompose the H₂ Hamiltonian into an appropriate operation
    /// that we can apply to qubits as we please.
    operation H2CompiledTrotterStep (idxBondLength : Int, trotterOrder : Int, trotterStepSize : Double, qubits : Qubit[], evGen : EvolutionGenerator) : Unit is Adj + Ctl {
        let simulationAlgorithm = TrotterSimulationAlgorithm(trotterStepSize, trotterOrder);
        simulationAlgorithm!(trotterStepSize, evGen, qubits);
    }

    // From a generator, converts it into a generator describing a compiled Hamiltonian
    operation ConvertGenerator(newMaxNum : Int, gen : EvolutionGenerator) : EvolutionGenerator {
        let (evolSet, genSys) = gen!;
        let (maxNum, genIndexGen) = genSys!;

        // constants
        // let newMaxNum = ????;

        // let's get each coefficient
        mutable coeffs = new Double[0];
        mutable totalCoeff = 0.0;
        for (index in 0..maxNum - 1) {
            // let's get the generator index
            let tempIndexSys = genIndexGen(index);
            let ((pauliTypes, coeff), targets) = tempIndexSys!;

            set totalCoeff += AbsD(coeff[0]);
            set coeffs += [coeff[0]];
        }

        // now, let's randomly select our newMaxNum types to apply
        mutable elementsToApply = new Int[0];
        for (idx in 0..newMaxNum - 1) {
            // We need to randomly select one from the coeffs
            mutable element = Random(Mapped(AbsD(_), coeffs));
            // Message($"{Random([1.0, 2.0, 2.0])}");
            // Message($"{Random([1.0, 2.0, 2.0])}");
            // Message($"{Random([1.0, 2.0, 2.0])}");

            // correct Q# bug
            if (element == maxNum) {
                set element = maxNum - 1;
            }


            Fact(element < maxNum, $"{element}, {maxNum} the max element");
            Fact(element >= 0, $"{element}, {maxNum} the min element");

            set elementsToApply += [element];
        }

        // Int -> GeneratorIndex
        let tau = totalCoeff / IntAsDouble(newMaxNum);
        let newGenIndex = ReturnAppropriateElement(_, elementsToApply, genIndexGen, coeffs, tau);
        let newGenSys = GeneratorSystem(newMaxNum, newGenIndex);
        return EvolutionGenerator(evolSet, newGenSys);
    }

    // given an index, looks up the element it should be mapped to, and plugs the value into the genIndex
    // startIdx = 0, mapArr = [1, 3, 2], genIndex 
    // returns genIndex(1)
    function ReturnAppropriateElement(startIdx : Int, mapArr : Int[], genIndex : (Int -> GeneratorIndex), coeffs : Double[], newCoeff : Double) : GeneratorIndex {
        let newIdx = mapArr[startIdx];
        
        // unpack the old idx
        let orig = genIndex(newIdx);
        let ((types, coeff), targets) = orig!;

        // correct for the coefficient if its negative
        let signedCoeff = newCoeff * IntAsDouble(SignD(coeffs[newIdx]));
        return GeneratorIndex((types, [signedCoeff]), targets);
    }
}


