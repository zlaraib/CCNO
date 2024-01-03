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
}}



	//=======//
	// Tests //
	//=======//
	stage('Basic Expectation Values'){ steps{
		sh 'julia tests/test_file.jl'
    } 
}
	stage('Particles evolution'){ steps{
		sh 'julia tests/par_evolution.jl'
		archiveArtifacts artifacts: 'misc/plots/evol/*/*/*.pdf'
    } 
}
	stage('Only Vacuum oscillations'){ steps{
		sh 'julia tests/main_vac_osc.jl'
		archiveArtifacts artifacts: 'misc/plots/vac_osc/*/*/*.pdf'
    }
} 
	stage('Rogerro(2021)_only_self_interactions'){ steps{
		sh 'julia tests/main_self_interaction_Rog.jl'
		archiveArtifacts artifacts: 'misc/plots/Rog_self_int/*/*/*.pdf'
    } 
}
	stage('Richers(2021)_Homogenous_FFI'){ steps{
		sh 'julia tests/Homogenous_FFI_Richers.jl'
		archiveArtifacts artifacts: 'misc/plots/FFI/*/*/*.pdf'
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
