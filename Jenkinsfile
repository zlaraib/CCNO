pipeline {
    triggers { pollSCM('') }  // Run tests whenever a new commit is detected.
    agent { dockerfile {args '--gpus all'}} // Use the Dockerfile defined in the root Flash-X directory
    stages {

        //=============================//
    	// Set up submodules and amrex //
        //=============================//
    	stage('Prerequisites'){ steps{
	    sh 'mpicc -v'
	    sh 'nvidia-smi'
	    sh 'nvcc -V'
	    sh 'git submodule update --init'
	    sh 'julia -v'

	    // use the Manifest.toml and Project.toml to install prerequisites
	    sh 'julia -e \'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()\''
}}



	//=======//
	// Tests //
	//=======//
	stage('Only Vacuum oscillations'){ steps{
		sh 'rm -rf test/datafiles; julia --project=. test/main_vac_osc.jl'
		archiveArtifacts artifacts: '*.pdf'
    }
} 
	stage('Richers(2021)_Bipolar_Oscillations'){ steps{
		sh 'rm -rf test/datafiles; julia --project=. test/Bipolar_Osc_Richers.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Richers(2021)_Homogenous_FFI'){ steps{
		sh 'rm -rf test/datafiles; julia --project=. test/Homogenous_FFI_MF_Richers.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Richers(2021)_Inhomogenous_FFI'){ steps{
		sh 'rm -rf test/datafiles; julia --project=. test/Inhomogenous_FFI_MF_Richers.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Roggero(2021)_only_self_interactions'){ steps{
		sh 'rm -rf test/datafiles; julia --project=. test/main_self_interaction_Rog.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Roggero(2021) full Hamiltonian'){ steps{
		sh 'rm -rf test/datafiles; julia --project=. test/main_Roggero.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
}// stages{

    post {
        always {
	    cleanWs(
	        cleanWhenNotBuilt: true,
		deleteDirs: true,
		disableDeferredWipeout: false,
		notFailBuild: true,
		patterns: [[pattern: 'amrex', type: 'EXCLUDE']] ) // allow amrex to be cached
	}
    }

} // pipeline{