#' Title
#'
#' @param theta
#' @param betX
#' @param pi1
#' @param sig1
#' @param w8s
#' @param m0
#' @param nX
#' @param bn
#' @param bins
#' @param M
#'
#' @importFrom stats fft sd
#'
#' @return
#' NOT EXPORTED @export
#' @keywords internal
#'
#' @examples
singleTrait_likelihood <- function(theta,betX,pi1,sig1,w8s,m0,nX,bn=2^7,bins=10,M=1e7){

  M = M #1e7
  piX = abs(theta[1]);
  h2X = abs(theta[2]);
  iX = abs(theta[3]);

  sigX = sqrt(h2X/(piX*M))
  # Number of genotyped SNPs
  m = length(betX)

  Rp = iX/nX
  bX = array(betX, c(1,m)) # reshape(betXY(:,1),[1,1,m]);

  # Define grid for FFT
  minX = mean(bX)-(5*sd(bX));
  maxX = mean(bX)+(5*sd(bX));
  dX = (maxX-minX)/(bn-1);
  minX = minX-dX/2;
  maxX = maxX+dX/2;

  bXi = ceiling((bX-minX)/dX);
  bXi[bXi<1] = 1;
  bXi[bXi>bn] = bn;


  if(piX > 0.2 || piX < 1e-6 || h2X > 1 || h2X < 1e-6 || iX > 1.5 || iX < 0.5 ){
    logL = 1e10
  }else{
    min_pi1 = min(pi1)-1e-10;
    max_pi1 = max(pi1)+1e-10;
    dp = (max_pi1-min_pi1)/bins;
    pc = min_pi1 + (dp * matrix(seq(0.5,(bins-0.5),1),ncol=1))
    pix = ceiling((pi1-min_pi1)/dp);
    min_sig1 = min(sig1)-1e-10;
    max_sig1 = max(sig1)+1e-10;
    ds = (max_sig1-min_sig1)/bins;
    sc = min_sig1 + (ds * matrix(seq(0.5,(bins-0.5),1),ncol=1))
    six = ceiling((sig1-min_sig1)/ds);
    cix = pix + bins*(six-1);
    uni_cix = sort(unique(cix))
    ucix = match(uni_cix, cix)
    ixMap = match(cix, uni_cix)
    Sig1 = sc[six[ucix]];
    Pi1 = pc[pix[ucix]];
    mm = length(Sig1);

    Ax = aperm( array(rep((1/sigX) / Sig1, bn), c(mm,bn)), c(2,1))
    Qx = aperm( array(rep(Pi1 * piX, bn), dim = c(mm, bn)),c(2,1))


    j = 0:(bn-1)
    vi = 2*pi*(j-bn/2)/(maxX-minX-dX);

    Rx = array(rep(vi, mm), c(bn,mm))

    Lx = -m0 * ( 1 - 1 / sqrt(1 + Rx^2/Ax^2))*Qx
    Le = -(1/2)*(Rp*Rx^2)

    # /!\ complex numbers here!
    mf_init = -2 * log(as.complex(-1)) * ( (minX+dX/2) / (maxX-minX-dX) )*j
    mf = array(rep(mf_init, mm), dim=c(bn, mm))

    phi = exp(Lx+Le+mf);

    # In R, not possible to chose dimension for FFT, so we need a loop to do it for all rho bins
    FFT=array(NA, dim=c(bn, mm))
    for(l in 1:mm){
      FFT[,l] = fft(phi[,l])
    }

    FFTmod_init = (1/(maxX-minX-dX))*(as.complex(-1))^(bn*((minX+dX/2)/(maxX-minX-dX)) + j)
    FFTmod = array(rep(FFTmod_init, mm), dim=c(bn,mm))

    FFT0 = Re(FFT*FFTmod)
    pfE = FFT0[cbind(t(bXi), ixMap)];
    length(which(pfE<0))
    # If some, remove them & update weights to keep the correct set of SNPs
    my_w8s = w8s[pfE>0]
    pfE=pfE[pfE>0]

    # We use m * mean(...) to account for SNPs that may have been excluded before
    logL = -m * mean(log(pfE*my_w8s))
  }
  return(logL)
}

