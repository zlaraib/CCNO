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
	stage('Particles evolution'){ steps{
		sh 'julia tests/par_evolution.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Only Vacuum oscillations'){ steps{
		sh 'julia tests/main_vac_osc.jl'
		archiveArtifacts artifacts: '*.pdf'
    }
} 
	stage('Rogerro(2021)_only_self_interactions'){ steps{
		sh 'julia tests/main_self_interaction_Rog.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Rogerro(2021) full Hamiltonian'){ steps{
		sh 'julia tests/main_Rogerro.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Rogerro(2021) Bipolar'){ steps{
		sh 'julia tests/main_Bipolar_Rog.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Rog_full H looped over N'){ steps{
		sh 'julia tests/main_Rog_N_loop.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('t_p vs symmetric delta_omega (Rog)'){ steps{
		sh 'julia tests/t_p_vs_sym_delta_w.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('t_p vs N (unsymmetric)(Rog)'){ steps{
		sh 'julia tests/t_p_vs_N_unsym.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('t_p vs N (symmetric)(Rog)'){ steps{
		sh 'julia tests/t_p_vs_N_sym.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Richers(2021)_Bipolar_Oscillations'){ steps{
		sh 'julia tests/Bipolar_Osc_Richers.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Richers(2021)_Homogenous_FFI'){ steps{
		sh 'julia tests/Homogenous_FFI_MF_Richers.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Homogenous_FFI_MB'){ steps{
		sh 'julia tests/Homo_FFI_MB.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Richers(2021)_Inhomogenous_FFI'){ steps{
		sh 'julia tests/Inhomogenous_FFI_MF_Richers.jl'
		archiveArtifacts artifacts: '*.pdf'
    } 
}
	stage('Inhomogenous_FFI_MB'){ steps{
		sh 'julia tests/Inhomo_FFI_MB.jl'
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