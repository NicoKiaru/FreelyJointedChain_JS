/*
 * FJC.java
 *
 * Created on 9 fevrier 2008, 01:02
 *
 * @author nicolas chiaruttini: nicolas.chiaruttini {at} gmail.com
 * translation from java applet to javascript started on 02/07/2012
 * GitHub Transfer 30/10/2017
 */

var simul_dt=0.01; // time step for integration

function FJC(){
	this.N=200; // monomer numbers
	this.d=0.02; // spacing between monomers
	this.masse=0.1; //weigth per monomer
	
	this.px=new Array(this.N); // X position array at t
	this.py=new Array(this.N); // Y position array at t
	this.pz=new Array(this.N); // Z position array at t

	this.bouge=new Array(this.N); // flags if monomer is allowed to move
	this.gamma=1; 	// link with temperature reservoir
	this.Theta=2;     // Temperature
	
	this.ErrTolerance=0.01;
	// precomputed numbers for speed purpose
	this.dcarre=this.d*this.d; 
	this.usm=1/this.masse; 
    	this.usdt=1/simul_dt
	this.dtc=simul_dt*simul_dt;

	this.ecart_type_force; //strength of brownian forces, depend on gamma and theta

	this.pxa=new Array(this.N); // X position array at t-dt
	this.pya=new Array(this.N); // Y position array at t-dt
	this.pza=new Array(this.N); // Z position array at t-dt
	
	this.pxu=new Array(this.N); // X position array at t, unconstrained then constrained
	this.pyu=new Array(this.N); // Y position array at t, unconstrained then constrained
	this.pzu=new Array(this.N); // Z position array at t, unconstrained then constrained	

	this.vx=new Array(this.N); // X speed at t
	this.vy=new Array(this.N); // Y speed at t
	this.vz=new Array(this.N); // Z speed at t

	this.fx=new Array(this.N); // X force at t
	this.fy=new Array(this.N); // Y force at t
	this.fz=new Array(this.N); // Z force at t
}

FJC.prototype.Init = function() {
	this.ecart_type_force=Math.sqrt(4.0*this.Theta*this.usm*this.gamma/simul_dt);
	for (var i=0; i<this.N; i++) {
		this.bouge[i]=true;
	        this.px[i]=0;
	        this.py[i]=this.d*i-this.d*(this.N/2.0);
	        this.pz[i]=0;
	        this.pxa[i]=0;
	        this.pya[i]=this.py[i];
	        this.pza[i]=0;
	        this.pxu[i]=0;
	        this.pyu[i]=this.py[i];
	        this.pzu[i]=0;
	        this.vx[i]=0;this.vy[i]=0;this.vz[i]=0;
	       this.fx[i]=0;this.fy[i]=0;this.fz[i]=0;
	}
}
    
FJC.prototype.Rand_Gauss = function() {
	return ((Math.random()+Math.random()+Math.random()+Math.random()-2)/4*this.ecart_type_force);
}
    
FJC.prototype.CalculForces = function(){
	  for (var i=0; i<this.N; i++) {
          	// random force -> brownian motion
          	this.fx[i]=this.masse*this.Rand_Gauss();
          	this.fy[i]=this.masse*this.Rand_Gauss();
          	this.fz[i]=this.masse*this.Rand_Gauss();
          	// drag force / dissipation
          	this.fx[i]+=this.masse*this.gamma*this.usdt*(this.pxa[i]-this.px[i]);
          	this.fy[i]+=this.masse*this.gamma*this.usdt*(this.pya[i]-this.py[i]);
          	this.fz[i]+=this.masse*this.gamma*this.usdt*(this.pza[i]-this.pz[i]);         
          	//this.fy[i]+=force_champ;
     	   }
}
    
FJC.prototype.CalculPositions = function() {
	var i,j;
	var errmax,lambda,erract;
	// Verlet integration
	// Computing non constrained new positions
	this.ecart_type_force = Math.sqrt(4.0*this.Theta*this.usm*this.gamma/simul_dt);
	//t=t+dt;

	for (var i=0; i<this.N; i++) {
        	if (this.bouge[i]==true) {
        		this.pxu[i]=2.0*this.px[i]-this.pxa[i]+this.usm*this.fx[i]*this.dtc;
        		this.pyu[i]=2.0*this.py[i]-this.pya[i]+this.usm*this.fy[i]*this.dtc;
	        	this.pzu[i]=2.0*this.pz[i]-this.pza[i]+this.usm*this.fz[i]*this.dtc;
		 }
	}

	// Constrains handling : 
	//    - fixed distance between monomers :
	//          SHAKE iterative method, 
	//    - polymer should not cross borders... (neglected here)
	errmax=1;
        while (errmax>this.ErrTolerance){
         	errmax=0;
       	 	// first - links
         	for (i=0; i<this.N-1; i++) { // for each link			
         		erract=(dist_carre(this.pxu,this.pyu,this.pzu,i,i+1)-this.dcarre); // current error
        		if ((erract*erract)>errmax) {errmax=erract*erract;}        
         		lambda=-.25* erract/prod_scal(this.px,this.py,this.pz,this.pxu,this.pyu,this.pzu,i,i+1);
        		// - pxu et pyu recomputation
       			if (this.bouge[i]){
       	  	              this.pxu[i]=this.pxu[i]+lambda*(this.px[i]-this.px[i+1]);
       	  	              this.pyu[i]=this.pyu[i]+lambda*(this.py[i]-this.py[i+1]);
       	          	      this.pzu[i]=this.pzu[i]+lambda*(this.pz[i]-this.pz[i+1]);
		  	}
                  	if (this.bouge[i+1]) {
       	                	this.pxu[i+1]=this.pxu[i+1]-lambda*(this.px[i]-this.px[i+1]);
       	                	this.pyu[i+1]=this.pyu[i+1]-lambda*(this.py[i]-this.py[i+1]);
       	                	this.pzu[i+1]=this.pzu[i+1]-lambda*(this.pz[i]-this.pz[i+1]);
       		  	}
        		// next walls (nothing here)
   		}
                // position updates : pxa<-px, et px<-pxu   
   		for (i=0; i<this.N; i++) {
       	   		if (this.bouge[i]) {
       				this.pxa[i]=this.px[i];this.pya[i]=this.py[i];this.pza[i]=this.pz[i];
       				this.px[i]=this.pxu[i];this.py[i]=this.pyu[i];this.pz[i]=this.pzu[i];
	   		}
   		}
	}
}


function prod_scal(x1,y1,z1,x2,y2,z2,i,j){
	return ((x1[j]-x1[i])*(x2[j]-x2[i])+(y1[j]-y1[i])*(y2[j]-y2[i])+(z1[j]-z1[i])*(z2[j]-z2[i]));
}

// Returns squared distance between monomer i and j, data are contained in x, y, z arrays
function dist_carre(x,y,z,i,j){
	return (x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]);
}

window.requestAnimFrame = (function(){
    return window.requestAnimationFrame       || // Standardized call
           window.webkitRequestAnimationFrame || // Chrome and Safari
           window.mozRequestAnimationFrame    || // Firefox
           window.oRequestAnimationFrame      || // Opera
           window.msRequestAnimationFrame     || // Internet Explorer
           function(callback){                   // Others...
               window.setTimeout(callback, 1000 / 60); // 60 fps
           };
})();

window.onload = function() {
    var canvas  = document.querySelector('#canvas'); // Animates in canvas part
    var context = canvas.getContext('2d');
	context.lineWidth = "5";
	context.strokeStyle = "black";
	context.lineJoin = 'round';
	context.lineCap = 'round';
    var fjc= new FJC();
    fjc.Init();
    var sc=250;
    function draw(){
        context.clearRect(0, 0, 500, 500);		
		fjc.CalculForces();
		fjc.CalculPositions();
		context.beginPath();
		var cx=0;
		var cy=0;
		for (i=1; i<fjc.N; i++) {
			cx=cx+fjc.px[i];
			cy=cy+fjc.py[i];
		}
		cx=cx/fjc.N;
		cy=cy/fjc.N;
		context.moveTo((fjc.px[0]-cx)*sc+150, (fjc.py[0]-cy)*sc+150);
		for (i=1; i<fjc.N; i++) {
			context.lineTo((fjc.px[i]-cx)*sc+150, (fjc.py[i]-cy)*sc+150);
		}    
		context.stroke();
        window.requestAnimFrame(function() { draw() });
    }
    draw(); // First call
};
