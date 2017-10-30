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
	this.N=50; // monomer numbers
	this.d=0.02; // spacing between monomers
	this.masse=0.1; //weigth per monomer
	
	this.px=new Array(this.N); // X position array
	this.py=new Array(this.N); // Y position array
	this.pz=new Array(this.N); // Z position array

	this.bouge=new Array(this.N); // flags if monomer is allowed to move
	this.gamma=1; 	// link with temperature reservoir
	this.Theta=2;     // TempÃ©rature
	
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
	
	this.pxu=new Array(this.N); // X position array at ?
	this.pyu=new Array(this.N); // Y position array at ?
	this.pzu=new Array(this.N); // Z position array at ? 	

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
          	// partie alÃ©atoire
          	this.fx[i]=this.masse*this.Rand_Gauss();
          	this.fy[i]=this.masse*this.Rand_Gauss();
          	this.fz[i]=this.masse*this.Rand_Gauss();
          	// partie dissipative
          	this.fx[i]+=this.masse*this.gamma*this.usdt*(this.pxa[i]-this.px[i]);
          	this.fy[i]+=this.masse*this.gamma*this.usdt*(this.pya[i]-this.py[i]);
          	this.fz[i]+=this.masse*this.gamma*this.usdt*(this.pza[i]-this.pz[i]);         
          	//this.fy[i]+=force_champ;

     	   }
}
    
FJC.prototype.CalculPositions = function() {
	var i,j;
	var errmax,lambda,erract;
	// IntÃ©gration par l'algorithme de Verlet
	// Calcul des nouvelles positions non contraintes
	this.ecart_type_force = Math.sqrt(4.0*this.Theta*this.usm*this.gamma/simul_dt);
	//t=t+dt;

	for (var i=0; i<this.N; i++) {
        	if (this.bouge[i]==true) {

        		this.pxu[i]=2.0*this.px[i]-this.pxa[i]+this.usm*this.fx[i]*this.dtc;
        		this.pyu[i]=2.0*this.py[i]-this.pya[i]+this.usm*this.fy[i]*this.dtc;
	        	this.pzu[i]=2.0*this.pz[i]-this.pza[i]+this.usm*this.fz[i]*this.dtc;
		 }
	}

	// Prise en compte des contraintes : 
	//    - la distance entre chaque monomÃ¨re reste constante :
	//          On utilise la mÃ©thode itÃ©rative SHAKE, 
	//    - le polymÃ¨re ne doit pas traverser les murs...
	errmax=1;
        while (errmax>this.ErrTolerance){
         	errmax=0;
       	 	// d'abord les liens
         	for (i=0; i<this.N-1; i++) { // Pour chaque lien
			
         		erract=(dist_carre(this.pxu,this.pyu,this.pzu,i,i+1)-this.dcarre);
        		if ((erract*erract)>errmax) {errmax=erract*erract;}        
         		lambda=-.25* erract/prod_scal(this.px,this.py,this.pz,this.pxu,this.pyu,this.pzu,i,i+1);
        		// - on recalcule pxu et pyu
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
        		  	// ensuite les murs...
   		}
                // on actualise pxa<-px, et px<-pxu   
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

    // Fonction qui renvoie la distance au carre entre deux points
    // i et j dont les coordonnees sont dans les tableaux x et y
function dist_carre(x,y,z,i,j){
	return (x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]);
}

window.requestAnimFrame = (function(){
    return window.requestAnimationFrame       || // La forme standardisÃ©e
           window.webkitRequestAnimationFrame || // Pour Chrome et Safari
           window.mozRequestAnimationFrame    || // Pour Firefox
           window.oRequestAnimationFrame      || // Pour Opera
           window.msRequestAnimationFrame     || // Pour Internet Explorer
           function(callback){                   // Pour les Ã©lÃ¨ves du dernier rang
               window.setTimeout(callback, 1000 / 60);
           };
})();

window.onload = function() {
    var canvas  = document.querySelector('#canvas');
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
    draw(); // premier appel
};
