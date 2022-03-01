#' Title
#'
#' @param theta
#' @param betXY
#' @param pi1
#' @param sig1
#' @param w8s
#' @param m0
#' @param pi_U
#' @param pi_X
#' @param pi_Y
#' @param i_X
#' @param i_Y
#' @param nX
#' @param nY
#' @param model
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
pairTrait_twoStep_likelihood <- function(theta,betXY,pi1,sig1,w8s,m0,pi_U=0.1,pi_X,pi_Y,i_X,i_Y,nX,nY,model="comp",bn=2^8,bins=15,M=1e7){

  M = M #1e7
  piX = pi_X
  piU = pi_U #piU is unidentifiable
  piY = pi_Y
  h2X = abs(theta[1]);
  h2Y = abs(theta[2]);
  iX = i_X;
  iY = i_Y;

  if(model=="comp"){
    tX = abs(theta[3]);
    tY = theta[4];
    axy = theta[5];
    ayx = theta[6];
    iXY = theta[7];
  }

  # if(model=="axy"){
  #   tX = abs(theta[3]);
  #   tY = theta[4];
  #   axy = 0;
  #   ayx = theta[5];
  #   iXY = theta[6];
  # }
  #
  # if(model=="ayx"){
  #   tX = abs(theta[3]);
  #   tY = theta[4];
  #   axy = theta[5];
  #   ayx = 0;
  #   iXY = theta[6];
  # }
  #
  # if(model=="tX"){
  #   tX = 0;
  #   tY = theta[3];
  #   axy = theta[4];
  #   ayx = theta[5];
  #   iXY = theta[6];
  # }
  #
  # if(model=="tY"){
  #   tX = abs(theta[3]);
  #   tY = 0;
  #   axy = theta[4];
  #   ayx = theta[5];
  #   iXY = theta[6];
  # }

  if(model=="U"){
    tX = 0;
    tY = 0;
    axy = theta[3];
    ayx = theta[4];
    iXY = theta[5];
  }


  sigX = sqrt(h2X/(piX*M))
  sigU = sqrt(1/(piU*M)); #h2U = 1
  sigY = sqrt(h2Y/(piY*M))

  m = nrow(betXY)

  Rp = matrix(c(iX/nX,iXY/sqrt(nX*nY),iXY/sqrt(nX*nY),iY/nY), nrow=2, ncol=2, byrow = T)
  bX = array(betXY[,1], c(1,1,m))
  bY = array(betXY[,2], c(1,1,m))

  minX = mean(bX)-(5*sd(bX));
  maxX = mean(bX)+(5*sd(bX));
  minY = mean(bY)-(5*sd(bY));
  maxY = mean(bY)+(5*sd(bY));

  dX = (maxX-minX)/(bn-1);
  dY = (maxY-minY)/(bn-1);

  minX = minX-dX/2;
  maxX = maxX+dX/2;
  minY = minY-dY/2;
  maxY = maxY+dY/2;

  bXi = ceiling((bX-minX)/dX);
  bYi = ceiling((bY-minY)/dY);

  bXi[bXi<1] = 1;
  bXi[bXi>bn] = bn;
  bYi[bYi<1] = 1;
  bYi[bYi>bn] = bn;

  if(max(c(piX,piU,piY)) > 0.2 || min(c(piX,piU,piY)) < 1e-6 || max(c((h2X+tX^2),(h2Y+tY^2)))>1
     || min(c((h2X+tX^2),(h2Y+tY^2))) < 1e-6|| abs(axy)>=1 || abs(ayx)>=1 || min(c(iX,iY))<=0.5
     || max(c(iX,iY))>=1.5 || abs(iXY)>1){
    logL = 1e10
    nmiss = 0
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

    sigK = aperm(array(rep(Sig1, bn*mm), c(mm,bn,bn)), c(2,3,1))
    piK = aperm(array(rep(Pi1, bn*mm), c(mm,bn,bn)), c(2,3,1))

    j = (0:(bn-1))
    vi = 2*pi*(j-bn/2)/(maxX-minX-dX);
    wj = 2*pi*(j-bn/2)/(maxY-minY-dY);
    Rx = array(rep(vi, bn*mm), dim = c(bn, bn, mm))
    Ry = aperm( array(rep(wj, bn*mm), dim = c(bn, bn, mm)),c(2,1,3))

    coX = (Rx+axy*Ry)*sigK;
    coY = (Ry+ayx*Rx)*sigK;
    coU = sigU * ((tX+ayx*tY)*Rx + (tY+axy*tX)*Ry) * sigK;

    Lx = -m0*piX*( 1 - 1/sqrt(1+sigX^2*coX^2) )*piK;
    Ly = -m0*piY*( 1 - 1/sqrt(1+sigY^2*coY^2) )*piK;
    Lu = -m0*piU*( 1 - 1/sqrt(1+coU^2) )*piK;
    Le = -(1/2)*( (Rp[1,1]*Rx^2+Rp[2,2]*Ry^2) + 2*Rp[1,2]*(Rx*Ry) );
    # /!\ complex numbers here!
    mf_init = -2*log(as.complex(-1)) *
      ( matrix( rep( ((minX+dX/2)/(maxX-minX-dX))*j, length(j) ), nrow = length(j) , byrow=T) +
          ( matrix( rep( ((minY+dY/2)/(maxY-minY-dY))%*%t(j), length(j) ),  nrow = length(j) ) ) )
    mf = array(rep(mf_init, mm), dim=c(bn, bn, mm))
    phi = exp(Lx+Ly+Lu+Le+mf);

    FFT=array(NA, dim=c(bn, bn, mm))
    for(l in 1:mm){
      FFT[,,l] = fft(phi[,,l])
    }

    FFTmod_init = (1/((maxX-minX-dX)*(maxY-minY-dY)))*(as.complex(-1))^(bn*((minX+dX/2)/(maxX-minX-dX) +
                                                                              (minY+dY/2)/(maxY-minY-dY)) + (outer(j,j,FUN="+")) )

    FFTmod = array(rep(FFTmod_init, mm), dim=c(bn,bn,mm));
    FFT0 = Re(FFT*FFTmod);
    ixF = cbind(bXi, bYi, ixMap);
    pfE = FFT0[ixF];
    selu = which(pfE>0);
    pfE = pfE[selu];
    my_w8s = w8s[selu]  #update weights vector if some SNPs are excluded
    nmiss = length(selu);

    lpfE = log(pfE);
    logL = -m * mean(log(pfE*my_w8s))
  }
  return(logL)
}
