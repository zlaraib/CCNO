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
		sh 'julia tests/exp_file.jl'
    } 
}
	stage('Time Evolution'){ steps{
		sh 'julia tests/time_evol.jl'
    }
} 
	stage('Rog_main_serial'){ steps{
		sh 'julia tests/Rog_tests/main_serial.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Rog_expect_val'){ steps{
		sh 'julia tests/Rog_tests/src/expect.jl'
    } 
}
	stage('Rog_create_gates'){ steps{
		sh 'julia tests/Rog_tests/src/gates_function.jl'
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
