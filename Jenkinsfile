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
	stage('Vacuum oscillation'){ steps{
		sh 'julia tests/main_vac_osc.jl'
		archiveArtifacts artifacts: '*.pdf'
    }
} 
	stage('Self interaction'){ steps{
		sh 'julia tests/main_self_interaction.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Rogerro(2021)_file'){ steps{
		sh 'julia tests/main_Rogerro.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Rog_particle_loop_test_file'){ steps{
		sh 'julia tests/Rog_particle_loop.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('t_p vs del_w (N=4) test_file'){ steps{
		sh 'julia tests/t_p_vs_delta_w_for_symmetric_w_a_and_w_b_using_N_4.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('t_p vs N(del_w=1, w_a=2, w_b=0) test_file'){ steps{
		sh 'julia tests/t_p_vs_N_for_w_a_2_w_b_0_using_del_w_1.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('t_p vs N(del_w=1, w_a=1, w_b=-1) test_file'){ steps{
		sh 'julia tests/t_p_vs_N_for_symmetric_w_a_and_w_b_using_del_w_1.jl'
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
