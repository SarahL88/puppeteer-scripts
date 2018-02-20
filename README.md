# Scripts for "PuppeteerLite" and "SD2 Puppeteer"
The scripts in this repository combine the SD2 demo workflow with the Living Computing Project workflow (aka PuppeteerLite) so that it works on the same UI

## Puppeteer UI
1) When a user uploads an SBOL file(or retrieves from SBH), the SD2 logic will be triggered     
Input: SBOL      
Output: Autoprotocol in json format   

2) When the user uploads a ZIP file, PuppeteerLite logic will be triggered     
*Important to note that PuppeteerLite was designed solely to meet the needs of the BUILD publication deadline*
Input: HeadtoHead2 archive    
Output: Tecan instructions in gwl format, human readable instructions in txt format    

##Flask Application
These scripts have been deployed on http://128.31.25.45:80 as a Flask RESTful service     

Deployment instructions https://www.digitalocean.com/community/tutorials/how-to-deploy-a-flask-application-on-an-ubuntu-vps

## Helpful commands
SSH command - `ssh ubuntu@128.31.25.45 -i bio-circuit`
Debug helpers - `journalctl | tail , check logs at `/var/log/apache2`

Testing Autoprotocol - `curl -H "Content-Type:application/json" -d @response.json -X POST "http://128.31.25.45:80/generate-autoprotocol" -o ap-response.json where response.json is the output from Puppeteer`      

Testing Puppeteer - curl -H "Content-Type: application/json" -d @request.json -X POST "http://128.31.25.45:8080/penguin-services-0.1.0/build" -o response.json`    
