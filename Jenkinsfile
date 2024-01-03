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
	stage('Rogerro(2021) Self interactions'){ steps{
		sh 'julia tests/main_self_interaction.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Rogerro(2021) full Hamiltonian'){ steps{
		sh 'julia tests/main_Rogerro.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Rog_full H looped over N'){ steps{
		sh 'julia tests/main_Rog_N_loop.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('t_p vs symmetric delta_omega'){ steps{
		sh 'julia tests/t_p_vs_sym_delta_w.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('t_p vs N (unsymmetric)'){ steps{
		sh 'julia tests/t_p_vs_N_unsym.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('t_p vs N (symmetric)'){ steps{
		sh 'julia tests/t_p_vs_N_sym.jl'
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
